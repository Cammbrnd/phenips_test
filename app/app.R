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
BRANCHE_GH     <- Sys.getenv("PHENIPS_BRANCHE", "main")   # "main" en production ; surchargeable pour tester une branche

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
                          "Potentiel thermique theorique : seules les mailles de 8 km qui croisent au moins",
                          "un peuplement de sapin / epicea (BD Foret IGN) sont affichees ; l'epicea present",
                          "dans un peuplement melange n'est pas detecte. Ce n'est pas un risque observe."),
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
                          "commence), pas forcement achevee. Seules les mailles de 8 km qui croisent au moins",
                          "un peuplement de sapin / epicea (BD Foret IGN) sont affichees."),
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
                          "pas une prediction validee : a confronter a un piegeage de controle.",
                          "Seules les mailles de 8 km qui croisent au moins un peuplement de pin",
                          "(BD Foret IGN) sont affichees."),
    couches = list(
      list(nom = "Generations", fichier = "steno_generations_{variante}_latest.geojson", champ_label = "libelle"),
      list(nom = "Essaimage",   fichier = "steno_essaimage_{variante}_latest.geojson",   champ_label = "libelle")
    )
  ),
  phenthau = list(
    label = "PHENTHAU (Thaumetopoea processionea, chene)",
    meta  = "phenthau_meta.json",
    aide  = paste("Stade atteint a la date des donnees ; L3 = debut du risque urticant (poils urticants).",
                  "Date d'eclosion = la plus precoce des 3 modeles (Custers, Wagenhoff, Meurisse).",
                  "Mortalite des L1 = famine entre l'eclosion et le debourrement du chene : elle n'est",
                  "determinee qu'une fois le debourrement atteint. La classe affichee est la borne haute de la",
                  "fourchette obtenue avec les 3 modeles d'eclosion (fourchette exacte au survol). Sources :",
                  "Halbig et al. 2024 (Forest Ecol. Manag. 552:121525) ; eclosion : Custers 2003,",
                  "Wagenhoff et al. 2014 (Eur. J. Entomol. 111:207-215), Meurisse et al. 2012",
                  "(Agric. Forest Entomol. 14:239-250)."),
    avertissement = paste("ATTENTION : modele PHENTHAUproc (Halbig et al. 2024) calibre sur Quercus robur en",
                          "Allemagne jusqu'a 800 m d'altitude. Hors de ce domaine (sud de la France, montagne),",
                          "les resultats sont une extrapolation, notamment la mortalite des L1 (eclosion tres",
                          "precoce par rapport au debourrement). Seules les mailles de 8 km qui croisent au",
                          "moins un peuplement de chenes decidus (BD Foret IGN) sont affichees."),
    couches = list(
      list(nom = "Stade (processionnaire du chene)", fichier = "phenthau_stades_latest.geojson",
           champ_label = "label"),
      list(nom = "Debut du risque urticant (L3)", fichier = "phenthau_l3_latest.geojson",
           champ_label = "date_l3"),
      list(nom = "Eclosion des L1", fichier = "phenthau_eclosion_latest.geojson",
           champ_label = "date_eclosion"),
      list(nom = "Debourrement du chene", fichier = "phenthau_debourrement_latest.geojson",
           champ_label = "date_debourrement"),
      list(nom = "Mortalite des L1 (famine)", fichier = "phenthau_mortalite_latest.geojson",
           champ_label = "label", champ_survol = "survol")
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


# =============================================================================
# HISTORIQUE 1960+ : tables annuelles par maille (data/historique/*.csv.gz)
# =============================================================================
url_historique <- function(nom) {
  glue("https://raw.githubusercontent.com/{UTILISATEUR_GH}/{REPO_GH}/{BRANCHE_GH}/data/historique/{nom}")
}

# prefixe_* : colonnes gen_<scenario> / doy_<scenario> ; col_* : colonnes uniques (CHAPY)
MODELES_H <- list(
  phenips = list(
    label = "PHENIPS-Clim (Ips typographus)", fichier = "phenips_annuel.csv.gz",
    scenarios = c("Ecorce exposee au sud" = "sunny", "Ecorce a l'ombre" = "shaded"),
    prefixe_gen = "gen", prefixe_doy = "doy", nom_envol = "premier essaimage",
    essence = "sapin / epicea",
    note = paste("Generations presentes au 31/12 : une generation est comptee des qu'elle a commence",
                 "(x.5 = couvee soeur).")),
  chapy = list(
    label = "CHAPY (Pityogenes chalcographus)", fichier = "chapy_annuel.csv.gz",
    scenarios = NULL, col_gen = "gen", col_doy = "doy_envol", nom_envol = "premier envol",
    essence = "sapin / epicea",
    note = paste("Generations presentes au 31/12 : une generation est comptee des qu'elle a commence",
                 "(x.5 = couvee soeur).")),
  steno = list(
    label = "STENO-DJ (Ips sexdentatus)", fichier = "steno_annuel.csv.gz",
    scenarios = c("Moyen (reference)" = "moyen", "Pire (developpement rapide)" = "pire",
                  "Moins pire (developpement lent)" = "moins_pire"),
    prefixe_gen = "gen", prefixe_doy = "doy", nom_envol = "premier envol",
    essence = "pin",
    note = "Generations achevees au 31/12 (maxima thermiques, temperature de l'air)."))

colonne_h <- function(m, scenario, quoi) {          # quoi : "gen" ou "doy"
  if (is.null(m$scenarios)) return(if (quoi == "gen") m$col_gen else m$col_doy)
  paste(if (quoi == "gen") m$prefixe_gen else m$prefixe_doy, scenario, sep = "_")
}

# Cache GLOBAL (partage entre sessions) : les tables ne changent qu'a la reconstruction.
cache_h <- new.env()

lire_table_h <- function(modele) {
  if (exists(modele, envir = cache_h, inherits = FALSE)) return(get(modele, envir = cache_h))
  d <- tryCatch({
    f <- tempfile(fileext = ".csv.gz")
    request(url_historique(MODELES_H[[modele]]$fichier)) |> req_perform(path = f)
    read.csv(gzfile(f))
  }, error = function(e) {
    message("Erreur lecture historique (", modele, ") : ", conditionMessage(e))
    NULL
  })
  if (is.null(d) || nrow(d) == 0) return(NULL)      # un echec n'est pas memorise
  d$cle <- d$x * 1e7 + d$y                           # identifiant de maille
  assign(modele, d, envir = cache_h)
  d
}

# Mailles de 8 km (carres L93 -> WGS84), triees par identifiant comme les matrices.
# Carres construits a la main : st_buffer(endCapStyle = "SQUARE") renvoie des cercles.
grille_h <- function(modele, d) {
  nom <- paste0("grille_", modele)
  if (exists(nom, envir = cache_h, inherits = FALSE)) return(get(nom, envir = cache_h))
  cel <- unique(d[, c("x", "y", "cle")])
  cel <- cel[order(cel$cle), ]
  carre <- function(x, y, h = 4000)
    sf::st_polygon(list(cbind(c(x - h, x + h, x + h, x - h, x - h), c(y - h, y - h, y + h, y + h, y - h))))
  g <- sf::st_sfc(mapply(carre, cel$x, cel$y, SIMPLIFY = FALSE), crs = 2154)
  g <- sf::st_transform(g, 4326)
  assign(nom, g, envir = cache_h)
  g
}

# Matrice mailles x annees d'une colonne
matrice_h <- function(d, col) {
  cles   <- sort(unique(d$cle))
  annees <- sort(unique(d$annee))
  M <- matrix(NA_real_, nrow = length(cles), ncol = length(annees), dimnames = list(NULL, annees))
  M[cbind(match(d$cle, cles), match(d$annee, annees))] <- d[[col]]
  list(M = M, cles = cles, annees = annees)
}

# Pente (par decennie) et p-value du test de correlation, maille par maille.
# Maille ignoree si moins de 20 annees renseignees ou valeur constante.
tendance_h <- function(M, annees) {
  r <- t(apply(M, 1, function(y) {
    ok <- !is.na(y); n <- sum(ok)
    if (n < 20 || length(unique(y[ok])) < 2) return(c(NA_real_, NA_real_))
    x <- annees[ok]; yy <- y[ok]
    xm <- mean(x)
    b  <- sum((x - xm) * (yy - mean(yy))) / sum((x - xm)^2)
    cc <- cor(x, yy)
    p  <- if (abs(cc) >= 1) 0 else 2 * pt(-abs(cc * sqrt((n - 2) / (1 - cc^2))), n - 2)
    c(b * 10, p)
  }))
  list(pente = r[, 1], p = r[, 2])
}

etiquette_doy <- function(type, cuts, p) format(as.Date(round(cuts) - 1, origin = "2001-01-01"), "%d/%m")


ui <- fluidPage(
  titlePanel("Phenologie forestiere - France metropolitaine"),
  tabsetPanel(
    id = "onglet",
    tabPanel(
      "Situation actuelle",
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
      leafletOutput("carte", height = "68vh")
    ),
    tabPanel(
      "Historique 1960-2025",
      fluidRow(
        column(3,
          selectInput("h_modele", "Modele",
                      choices = setNames(names(MODELES_H), sapply(MODELES_H, `[[`, "label"))),
          uiOutput("ui_h_scenario"),
          radioButtons("h_variable", "Variable",
                       choices = c("Nombre de generations" = "gen", "Date du premier envol" = "doy")),
          radioButtons("h_vue", "Affichage",
                       choices = c("Une annee" = "annee", "Moyenne sur une periode" = "moyenne",
                                   "Tendance (par decennie)" = "tendance")),
          conditionalPanel("input.h_vue == 'annee'",
            sliderInput("h_annee", "Annee", min = 1960, max = 2025, value = 2025, step = 1, sep = "",
                        animate = animationOptions(interval = 1000))),
          conditionalPanel("input.h_vue != 'annee'",
            sliderInput("h_periode", "Periode", min = 1960, max = 2025, value = c(1960, 2025),
                        step = 1, sep = "")),
          conditionalPanel("input.h_vue == 'tendance'",
            checkboxInput("h_signif", "Griser les tendances non significatives (p > 0,05)", value = TRUE)),
          uiOutput("note_h")
        ),
        column(9,
          uiOutput("bandeau_h"),
          leafletOutput("carte_h", height = "52vh"),
          plotOutput("serie_h", height = "230px")
        )
      )
    )
  ),
  tags$p(
    style = "font-size: 11px; color: #6c757d; margin: 6px 0 12px 0;",
    paste("Source : Meteo-France (donnees SAFRAN), via GeoSAS INRAE ; essences : IGN, BD Foret V2.",
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
    if (length(m$variantes) == 0) return(if (is.null(m$aide)) NULL else helpText(m$aide))   # PHENTHAU : aide sans scenario
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
      addProviderTiles("OpenStreetMap.France") |>
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
      # texte du survol distinct de la legende (optionnel, ex. fourchette de mortalite PHENTHAU)
      champ_survol <- premier_champ(g, cc$champ_survol)
      survol <- if (is.na(champ_survol)) labels else as.character(g[[champ_survol]])

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
                             weight = 0.5, fillOpacity = 0.65, group = cc$nom, label = survol)

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

  # ===========================================================================
  # ONGLET HISTORIQUE
  # ===========================================================================
  m_h <- reactive({ req(input$h_modele); MODELES_H[[input$h_modele]] })

  output$ui_h_scenario <- renderUI({
    m <- m_h()
    if (is.null(m$scenarios)) return(NULL)
    radioButtons("h_scenario", "Scenario", choices = m$scenarios, selected = unname(m$scenarios[1]))
  })

  output$note_h <- renderUI({
    m <- m_h()
    helpText(paste0("Mailles avec ", m$essence, " (BD Foret IGN). ", m$note),
             " Chaque annee est calculee de janvier a decembre, sans lien avec l'annee precedente.")
  })

  # Chargement a la demande : seulement quand l'onglet est ouvert
  table_h <- reactive({
    req(input$onglet == "Historique 1960-2025")
    lire_table_h(req(input$h_modele))
  })

  # Colonne choisie ; le scenario n'est lu que s'il appartient au modele courant
  col_h <- reactive({
    m <- m_h()
    sc <- if (is.null(m$scenarios)) NULL else
            if (isTRUE(input$h_scenario %in% m$scenarios)) input$h_scenario else unname(m$scenarios[1])
    colonne_h(m, sc, req(input$h_variable))
  })

  mat_h <- reactive({
    d <- table_h()
    shiny::validate(shiny::need(!is.null(d), "Tables historiques indisponibles (pas encore calculees ?)."))
    col <- col_h()
    shiny::validate(shiny::need(col %in% names(d), "Colonne absente de la table."))
    matrice_h(d, col)
  })

  # Les curseurs suivent les annees reellement presentes dans la table
  observeEvent(table_h(), {
    an <- range(table_h()$annee)
    updateSliderInput(session, "h_annee", min = an[1], max = an[2], value = an[2])
    updateSliderInput(session, "h_periode", min = an[1], max = an[2], value = an)
  })

  cle_choisie <- reactiveVal(NULL)
  observeEvent(list(input$h_modele, input$h_scenario, input$h_variable), cle_choisie(NULL), ignoreInit = TRUE)
  observeEvent(input$carte_h_shape_click, cle_choisie(as.numeric(input$carte_h_shape_click$id)))

  annee_h <- debounce(reactive(input$h_annee), 250)

  # Valeurs a cartographier : un vecteur par maille + palette + legende
  valeurs_h <- reactive({
    mt <- mat_h(); M <- mt$M
    vue <- req(input$h_vue); variable <- req(input$h_variable)
    domaine <- if (variable == "doy") quantile(M, c(0.01, 0.99), na.rm = TRUE) else range(M, na.rm = TRUE)
    domaine <- as.numeric(domaine)
    if (diff(domaine) == 0) domaine <- domaine + c(-0.5, 0.5)

    if (vue == "annee") {
      an <- req(annee_h())
      shiny::validate(shiny::need(as.character(an) %in% colnames(M), "Annee absente de la table."))
      v <- M[, as.character(an)]
      titre <- as.character(an)
    } else {
      per <- req(input$h_periode)
      cols <- which(mt$annees >= per[1] & mt$annees <= per[2])
      shiny::validate(shiny::need(length(cols) >= 2, "Choisir une periode d'au moins 2 annees."))
      if (vue == "moyenne") {
        v <- rowMeans(M[, cols, drop = FALSE], na.rm = TRUE); v[is.nan(v)] <- NA
        titre <- sprintf("moyenne %d-%d", per[1], per[2])
      } else {
        tr <- tendance_h(M[, cols, drop = FALSE], mt$annees[cols])
        v <- tr$pente
        titre <- sprintf("tendance %d-%d", per[1], per[2])
      }
    }

    if (vue == "tendance") {
      lim <- as.numeric(quantile(abs(tr$pente), 0.98, na.rm = TRUE)); if (!is.finite(lim) || lim == 0) lim <- 1
      pal <- colorNumeric("RdBu", domain = c(-lim, lim), na.color = "#cccccc",
                          reverse = (variable == "gen"))      # chaud = plus de generations / envol plus precoce
      unite <- if (variable == "gen") "generations par decennie" else "jours par decennie"
      list(v = pmax(pmin(v, lim), -lim), brut = v, p = tr$p, pal = pal, dom = c(-lim, lim),
           titre = titre, legende = unite, fmt = labelFormat(digits = 2), tendance = TRUE)
    } else {
      pal <- colorNumeric("YlOrRd", domain = domaine, na.color = "#cccccc", reverse = (variable == "doy"))
      list(v = pmax(pmin(v, domaine[2]), domaine[1]), brut = v, pal = pal, dom = domaine, titre = titre,
           legende = if (variable == "gen") "Generations" else "Premier envol (jj/mm)",
           fmt = if (variable == "gen") labelFormat(digits = 1) else etiquette_doy, tendance = FALSE)
    }
  })

  output$bandeau_h <- renderUI({
    vh <- tryCatch(valeurs_h(), error = function(e) NULL)
    if (is.null(vh)) return(NULL)
    div(style = "font-size: 13px; color: #495057; margin: 8px 0;",
        strong(paste0(m_h()$label, " - ", if (input$h_variable == "gen") "nombre de generations" else
                        paste0("date du ", m_h()$nom_envol), " - ", vh$titre)))
  })

  output$carte_h <- renderLeaflet({
    leaflet() |>
      addProviderTiles("OpenStreetMap.France") |>
      setView(lng = 2.5, lat = 46.5, zoom = 6)
  })

  # La carte n'existe qu'a l'ouverture de l'onglet : on attend qu'elle ait envoye ses limites
  carte_h_prete <- reactiveVal(FALSE)
  observeEvent(input$carte_h_bounds, carte_h_prete(TRUE), once = TRUE)

  observe({
    req(carte_h_prete())
    vh <- valeurs_h(); mt <- mat_h()
    g  <- grille_h(req(input$h_modele), table_h())
    v  <- vh$v
    couleur <- vh$pal(v)
    opac <- rep(0.75, length(v))
    if (isTRUE(vh$tendance) && isTRUE(input$h_signif)) {
      ns <- is.na(vh$p) | vh$p > 0.05
      couleur[ns] <- "#d9d9d9"; opac[ns] <- 0.45
    }
    texte <- if (isTRUE(vh$tendance)) {
      sprintf("%s %+.2f (p = %s)", vh$legende, vh$brut, ifelse(is.na(vh$p), "-", formatC(vh$p, digits = 2, format = "g")))
    } else if (input$h_variable == "doy") {
      ifelse(is.na(vh$brut), "pas d'envol", format(as.Date(round(vh$brut) - 1, origin = "2001-01-01"), "%d/%m"))
    } else {
      ifelse(is.na(vh$brut), "-", formatC(vh$brut, digits = 1, format = "f"))
    }
    proxy <- leafletProxy("carte_h") |> clearShapes() |> clearControls()
    proxy <- addPolygons(proxy, data = sf::st_sf(geometry = g), layerId = as.character(mt$cles),
                         fillColor = couleur, fillOpacity = opac, color = "#333333", weight = 0.3,
                         label = texte)
    addLegend(proxy, position = "bottomright", pal = vh$pal, values = vh$dom,
              title = vh$legende, labFormat = vh$fmt, opacity = 0.9)
  })

  output$serie_h <- renderPlot({
    mt <- mat_h(); M <- mt$M; an <- mt$annees
    moy <- colMeans(M, na.rm = TRUE)
    unite <- if (input$h_variable == "gen") "generations" else "jour de l'annee (1er envol)"
    i <- if (is.null(cle_choisie())) NA_integer_ else match(cle_choisie(), mt$cles)
    y <- if (is.na(i)) NULL else M[i, ]
    par(mar = c(3, 4.2, 2.2, 0.8))
    plot(an, moy, type = "l", lwd = 2, col = "grey55", xlab = "", ylab = unite,
         ylim = range(c(moy, y), na.rm = TRUE),
         main = if (is.null(y)) "Moyenne de toutes les mailles (cliquer sur une maille pour sa serie)"
                else "Maille choisie (rouge) et moyenne de toutes les mailles (gris)")
    ok <- is.finite(moy)
    if (sum(ok) > 3) abline(lm(moy[ok] ~ an[ok]), lty = 2, col = "grey30")
    if (!is.null(y)) {
      lines(an, y, col = "#C0392B", lwd = 2); points(an, y, col = "#C0392B", pch = 16, cex = 0.6)
      oky <- is.finite(y)
      if (sum(oky) > 3) abline(lm(y[oky] ~ an[oky]), lty = 2, col = "#C0392B")
    }
  })

}

shinyApp(ui, server)
