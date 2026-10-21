# =============================================================================
# App Shiny multi-modeles — visualisation seule
# -----------------------------------------------------------------------------
# Ne fait AUCUN calcul : lit les GeoJSON publies par GitHub Actions et les
# affiche sur une carte leaflet, avec un selecteur pour changer de modele.
#
# Chaque modele exporte deja son propre champ "couleur" (calcule cote R au
# moment de l'export) : l'app se contente de l'utiliser tel quel. Exception :
# les couches de dates (champ date_iso : onset PHENIPS, envol CHAPY) sont
# coloriees ici par un degrade a domaine FIXE.
#
# v2 :
#   [P7]  PHENIPS : selecteur de scenario d'exposition (sud / ombre)
#   [P9]  degrade de dates a domaine fixe (15/02 -> 31/08) : la meme date a
#         toujours la meme couleur, d'une semaine a l'autre
#   [P11] legende pour toutes les couches a classes
#   [P12] seule la premiere couche est visible par defaut
#   [P13] date des donnees + alerte de fraicheur
#   [P16] mention de la source (Meteo-France) et des modeles
# =============================================================================

library(shiny)
library(leaflet)
library(sf)
library(glue)
library(httr2)

UTILISATEUR_GH <- "Cammbrnd"
REPO_GH        <- "phenips_test"
BRANCHE_GH     <- "main"

# Au-dela de ce nombre de jours, la date des donnees est signalee en alerte.
# SAFRAN est mis a jour par lots : 5 a 26 jours de decalage constates.
SEUIL_ALERTE_JOURS <- 21L

url_geojson <- function(nom_fichier) {
  glue("https://raw.githubusercontent.com/{UTILISATEUR_GH}/{REPO_GH}/{BRANCHE_GH}/data/{nom_fichier}")
}

# --- Definition des modeles ---------------------------------------------------
# `couches` : nom affiche (= groupe leaflet), fichier "latest", champ(s) du
# libelle (le premier present est utilise pour l'infobulle ET la legende).
# {variante} dans un nom de fichier est remplace par le scenario choisi.
MODELES <- list(
  phenips = list(
    label     = "PHENIPS-Clim (Ips typographus)",
    variantes = c("Ecorce exposee au sud" = "sunny", "Ecorce a l'ombre" = "shaded"),
    meta      = "phenips_meta.json",
    aide      = paste("Aucun des deux scenarios n'est systematiquement le plus defavorable :",
                      "une ecorce exposee au sud depasse souvent 38.9 C en ete, seuil au-dela",
                      "duquel le developpement s'arrete (PHENIPS). En annee chaude, le scenario",
                      "a l'ombre peut donc donner autant, voire plus, de generations."),
    avertissement = paste("Une generation est comptee des qu'elle est presente (developpement",
                          "commence), pas forcement achevee : l'indice continu donne l'avancement.",
                          "Potentiel thermique theorique : la carte ne tient pas compte de la presence",
                          "d'epicea et n'indique pas un risque la ou il n'y a pas de population."),
    couches = list(
      list(nom = "Generations", fichier = "phenips_generations_{variante}_latest.geojson",
           champ_label = "label"),
      list(nom = "Indice continu (generations cumulees)",
           fichier = "phenips_indice_{variante}_latest.geojson", champ_label = "label"),
      list(nom = "Onset (1er essaimage)", fichier = "phenips_onset_{variante}_latest.geojson",
           champ_label = "date_onset")
    )
  ),
  chapy = list(
    label = "CHAPY (Pityogenes chalcographus)",
    avertissement = paste("Une generation est comptee des qu'elle est presente (developpement",
                          "commence), pas forcement achevee."),
    couches = list(
      list(nom = "Generations",   fichier = "chapy_generations_latest.geojson", champ_label = "label"),
      list(nom = "Premier envol", fichier = "chapy_onset_latest.geojson",
           champ_label = c("date_envol", "label"))
    )
  ),
  steno = list(
    label     = "STENO-DJ (Ips sexdentatus)",
    # Marges d'erreur de la these de Pineau (seuil 10.9 +/- 0.2 C, 506 +/- 7 DJ).
    # Le scenario "moyen" garde les noms de fichiers historiques (sans suffixe).
    variantes = c("Moyen (reference)" = "moyen", "Pire (developpement rapide)" = "pire",
                  "Moins pire (developpement lent)" = "moins_pire"),
    sans_suffixe = "moyen",
    aide      = paste("Les trois scenarios reprennent les marges d'erreur publiees (seuil",
                      "10.9 +/- 0.2 C ; 506 +/- 7 DJ par generation). Elles sont etroites :",
                      "elles ne couvrent ni la temperature sous ecorce, ni la maturation,",
                      "ni la photoperiode, ni la variabilite d'une annee a l'autre."),
    avertissement = paste("Generations = MAXIMA THERMIQUES (temperature de l'air, hotes toujours",
                          "disponibles, pas d'arret photoperiodique). Indicateur d'aide a la decision,",
                          "pas une prediction validee : a confronter a un piegeage de controle."),
    couches = list(
      list(nom = "Generations", fichier = "steno_generations_{variante}_latest.geojson", champ_label = "libelle"),
      list(nom = "Essaimage",   fichier = "steno_essaimage_{variante}_latest.geojson",   champ_label = "libelle")
    )
  ),
  monochamus = list(
    label = "Monochamus galloprovincialis (bivolt. A + B)",
    couches = list(
      list(nom = "Stade + bivoltisme A", fichier = "mono_stade_latest.geojson",
           champ_label = "label"),
      list(nom = "Bivoltisme B (gen2)",  fichier = "mono_bivoltB_gen2_latest.geojson",
           champ_label = "label")
    )
  )
)

TOUS_GROUPES  <- unique(unlist(lapply(MODELES, function(m) vapply(m$couches, `[[`, "", "nom"))))
MAX_COUCHES   <- max(vapply(MODELES, function(m) length(m$couches), 1L))
IDS_LEGENDES  <- paste0("legende_", seq_len(MAX_COUCHES))

# Champs servant a ordonner une legende, par ordre de preference
CHAMPS_ORDRE <- c("ordre", "borne_inf", "nb_generations", "stade", "statut_gen2",
                  "doy_debut", "generation")

# Etiquettes de la legende du degrade de dates : jours depuis 1970 -> jj/mm
# (labelFormat() de leaflet n'a pas de parametre date, d'ou cette fonction)
etiquette_date <- function(type, cuts, p) {
  format(as.Date(cuts, origin = "1970-01-01"), "%d/%m")
}

premier_champ <- function(g, candidats) {
  c_ok <- intersect(candidats, names(g))
  if (length(c_ok) == 0) NA_character_ else c_ok[1]
}

date_depuis_couche <- function(g) {
  for (champ in c("date_calcul", "date_situation", "date_mise_a_jour")) {
    if (champ %in% names(g)) {
      d <- suppressWarnings(as.Date(as.character(g[[champ]][1])))
      if (!is.na(d)) return(d)
    }
  }
  as.Date(NA)
}


ui <- fluidPage(
  titlePanel("Phenologie forestiere - France metropolitaine"),
  fluidRow(
    column(3,
      selectInput("modele", "Modele",
                  choices  = setNames(names(MODELES), sapply(MODELES, `[[`, "label")),
                  selected = "phenips"),
      uiOutput("ui_variante")
    ),
    column(9,
      uiOutput("bandeau"),
      uiOutput("avertissement")
    )
  ),
  leafletOutput("carte", height = "68vh"),
  tags$p(
    style = "font-size: 11px; color: #6c757d; margin: 6px 0 12px 0;",
    paste("Source : Meteo-France (donnees SAFRAN), via GeoSAS INRAE.",
          "Modeles : PHENIPS-Clim (package barrks, Jentschke 2025 ; PHENIPS : Baier et al. 2007),",
          "CHAPY (Ogris et al. 2020) ; STENO-DJ (Pineau et al. 2016) et Monochamus : modeles internes fondes sur la litterature.",
          "Les cartes representent un potentiel de developpement thermique theorique,",
          "pas une observation de terrain.")
  )
)

server <- function(input, output, session) {

  # --- Lecture / cache (par session) ------------------------------------------
  cache <- new.env()

  telecharger <- function(url, lecteur) {
    tryCatch({
      f <- tempfile()
      request(url) |> req_perform(path = f)
      lecteur(f)
    }, error = function(e) {
      message("Erreur lecture (", url, "): ", conditionMessage(e))
      NULL
    })
  }

  charger <- function(nom_fichier) {
    if (!exists(nom_fichier, envir = cache, inherits = FALSE)) {
      assign(nom_fichier,
             list(valeur = telecharger(url_geojson(nom_fichier),
                                       function(f) sf::st_read(f, quiet = TRUE))),
             envir = cache)
    }
    get(nom_fichier, envir = cache, inherits = FALSE)$valeur
  }

  charger_meta <- function(nom_fichier) {
    if (!exists(nom_fichier, envir = cache, inherits = FALSE)) {
      assign(nom_fichier,
             list(valeur = telecharger(url_geojson(nom_fichier),
                                       function(f) jsonlite::fromJSON(f))),
             envir = cache)
    }
    get(nom_fichier, envir = cache, inherits = FALSE)$valeur
  }

  # {variante} -> scenario choisi. Pour le scenario "sans suffixe" (ex. STENO moyen),
  # le "_{variante}" disparait en entier : steno_generations_latest.geojson.
  nom_fichier_couche <- function(cc, variante, sans_suffixe = NULL) {
    f <- cc$fichier
    if (is.null(variante)) return(f)
    if (!is.null(sans_suffixe) && identical(variante, sans_suffixe)) {
      return(gsub("_{variante}", "", f, fixed = TRUE))
    }
    gsub("{variante}", variante, f, fixed = TRUE)
  }

  # --- Selection courante (modele + scenario) -----------------------------------
  selection <- reactive({
    m <- MODELES[[input$modele]]
    v <- NULL
    if (length(m$variantes) > 0) {
      v <- if (isTRUE(input$variante %in% m$variantes)) input$variante else unname(m$variantes[1])
    }
    list(m = m, v = v)
  })

  output$ui_variante <- renderUI({
    m <- MODELES[[input$modele]]
    if (length(m$variantes) == 0) return(NULL)
    tagList(
      radioButtons("variante", "Scenario d'exposition",
                   choices = m$variantes, selected = unname(m$variantes[1])),
      helpText(m$aide)
    )
  })

  output$avertissement <- renderUI({
    m <- MODELES[[input$modele]]
    if (is.null(m$avertissement)) return(NULL)
    div(style = "font-size: 12px; color: #6c757d; margin: 4px 0 8px 0;", m$avertissement)
  })

  # --- [P13] Date des donnees et alerte de fraicheur --------------------------------
  date_donnees <- function(m, v) {
    if (!is.null(m$meta)) {
      meta <- charger_meta(m$meta)
      if (!is.null(meta$date_donnees)) {
        d <- suppressWarnings(as.Date(meta$date_donnees))
        if (!is.na(d)) return(d)
      }
    }
    for (cc in m$couches) {
      g <- charger(nom_fichier_couche(cc, v, m$sans_suffixe))
      if (!is.null(g)) {
        d <- date_depuis_couche(g)
        if (!is.na(d)) return(d)
      }
    }
    as.Date(NA)
  }

  output$bandeau <- renderUI({
    sel <- selection()
    d   <- date_donnees(sel$m, sel$v)
    base_style <- "padding: 6px 10px; border-radius: 4px; font-size: 13px; margin-top: 25px;"
    if (is.na(d)) {
      return(div(style = paste0(base_style, "background:#f8d7da; border:1px solid #dc3545; color:#58151c;"),
                 strong("Date des donnees inconnue"), " : fichiers indisponibles."))
    }
    age    <- max(0L, as.integer(Sys.Date() - d))
    alerte <- age > SEUIL_ALERTE_JOURS
    style  <- if (alerte) "background:#fff3cd; border:1px solid #e0a800; color:#664d03;"
              else        "background:#f1f3f5; border:1px solid #dee2e6; color:#495057;"
    div(style = paste0(base_style, style),
        strong(paste0("Donnees au ", format(d, "%d/%m/%Y"))),
        paste0(" (", age, " jours de decalage)"),
        if (alerte) " - donnees anciennes : SAFRAN est mis a jour par lots, le decalage peut atteindre plusieurs semaines.")
  })

  # --- Carte ----------------------------------------------------------------------
  output$carte <- renderLeaflet({
    leaflet() |>
      addProviderTiles("CartoDB.Positron") |>
      setView(lng = 2.5, lat = 46.5, zoom = 6)
  })

  observe({
    sel <- selection()
    m <- sel$m; v <- sel$v

    proxy <- leafletProxy("carte")
    for (id in IDS_LEGENDES) proxy <- removeControl(proxy, id)
    proxy <- clearGroup(proxy, TOUS_GROUPES)

    groupes <- character(0)
    for (k in seq_along(m$couches)) {
      cc <- m$couches[[k]]
      g  <- charger(nom_fichier_couche(cc, v, m$sans_suffixe))
      if (is.null(g) || nrow(g) == 0) next        # couche absente : rien d'ancien n'est conserve

      champ_label <- premier_champ(g, cc$champ_label)
      labels <- if (is.na(champ_label)) NULL else as.character(g[[champ_label]])

      if ("date_iso" %in% names(g)) {
        # [P9] degrade a domaine fixe 15/02 -> 31/08 de l'annee des donnees
        dates <- as.Date(g$date_iso)
        an <- suppressWarnings(as.integer(format(max(dates, na.rm = TRUE), "%Y")))
        if (is.na(an)) an <- as.integer(format(Sys.Date(), "%Y"))
        d_min <- as.numeric(as.Date(sprintf("%d-02-15", an)))
        d_max <- as.numeric(as.Date(sprintf("%d-08-31", an)))
        pal   <- colorNumeric("YlOrRd", domain = c(d_min, d_max), na.color = "#cccccc")
        valeurs <- pmin(pmax(as.numeric(dates), d_min), d_max)
        proxy <- addPolygons(proxy, data = g, fillColor = pal(valeurs), color = "#333333",
                             weight = 0.5, fillOpacity = 0.75, group = cc$nom, label = labels)
        proxy <- addLegend(proxy, position = "bottomright", pal = pal, values = c(d_min, d_max),
                           title = "Date (jj/mm)", labFormat = etiquette_date,
                           layerId = IDS_LEGENDES[k], group = cc$nom)
      } else {
        couleurs <- if ("couleur" %in% names(g)) g$couleur else "#31688E"
        proxy <- addPolygons(proxy, data = g, fillColor = couleurs, color = "#333333",
                             weight = 0.5, fillOpacity = 0.65, group = cc$nom, label = labels)

        # [P11] legende des classes : une ligne par couple (libelle, couleur)
        if ("couleur" %in% names(g) && !is.na(champ_label)) {
          d <- sf::st_drop_geometry(g)
          cle <- premier_champ(d, CHAMPS_ORDRE)
          if (!is.na(cle)) d <- d[order(d[[cle]]), ]
          leg <- unique(data.frame(label = as.character(d[[champ_label]]),
                                   couleur = as.character(d$couleur),
                                   stringsAsFactors = FALSE))
          proxy <- addLegend(proxy, position = "bottomleft", colors = leg$couleur,
                             labels = leg$label, title = cc$nom, opacity = 0.9,
                             layerId = IDS_LEGENDES[k], group = cc$nom)
        }
      }
      groupes <- c(groupes, cc$nom)
    }

    if (length(groupes) > 0) {
      proxy <- addLayersControl(proxy, overlayGroups = groupes,
                                options = layersControlOptions(collapsed = FALSE))
      # [P12] une seule couche visible au depart : les autres se superposeraient.
      # showGroup() explicite : un nom de groupe commun a deux modeles
      # ("Generations") resterait sinon masque si l'utilisateur l'avait coche off.
      proxy <- showGroup(proxy, groupes[1])
      if (length(groupes) > 1) proxy <- hideGroup(proxy, groupes[-1])
    }
  })
}

shinyApp(ui, server)
