# =============================================================================
# App Shiny multi-modeles — visualisation seule
# -----------------------------------------------------------------------------
# Ne fait AUCUN calcul : lit les GeoJSON publies par GitHub Actions et les
# affiche sur une carte leaflet, avec un selecteur pour changer de modele.
#
# Chaque modele exporte deja son propre champ "couleur" (calcule cote R au
# moment de l'export, dans xxx_functions.R) : l'app se contente de l'utiliser
# tel quel, sans reconstruire de palette. Ca garantit que la carte QGIS, la
# carte Shiny et n'importe quel autre client qui lit ces GeoJSON affichent
# exactement les memes couleurs.
#
# Exception : l'onset PHENIPS (date continue) n'a pas de champ "couleur"
# pre-calcule (une palette discrete ne rendrait pas justice a un degrade
# continu) -> degrade colorNumeric() calcule ici, cote client.
# =============================================================================

library(shiny)
library(leaflet)
library(sf)
library(glue)
library(httr2)

UTILISATEUR_GH <- "Cammbrnd"
REPO_GH        <- "phenips_test"
BRANCHE_GH     <- "main"

url_geojson <- function(nom_fichier) {
  glue("https://raw.githubusercontent.com/{UTILISATEUR_GH}/{REPO_GH}/{BRANCHE_GH}/data/{nom_fichier}")
}

# --- Definition des 4 modeles : nom affiche, fichiers "latest", champ label ---
MODELES <- list(
  phenips = list(
    label        = "PHENIPS-Clim (Ips typographus)",
    fichier_1    = "phenips_generations_latest.geojson",
    fichier_2    = "phenips_onset_latest.geojson",
    nom_couche_1 = "Generations",
    nom_couche_2 = "Onset (1er essaimage)",
    champ_label_1 = "label",
    champ_label_2 = "date_onset"
  ),
  chapy = list(
    label        = "CHAPY (Pityogenes chalcographus)",
    fichier_1    = "chapy_generations_latest.geojson",
    fichier_2    = "chapy_envol_latest.geojson",
    nom_couche_1 = "Generations",
    nom_couche_2 = "Envol",
    champ_label_1 = "label",
    champ_label_2 = "label"
  ),
  steno = list(
    label        = "STENO-DJ (Ips sexdentatus)",
    fichier_1    = "steno_generations_latest.geojson",
    fichier_2    = "steno_essaimage_latest.geojson",
    nom_couche_1 = "Generations",
    nom_couche_2 = "Essaimage",
    champ_label_1 = "libelle",
    champ_label_2 = "libelle"
  ),
  monochamus = list(
    label        = "Monochamus galloprovincialis (bivolt. A + B)",
    fichier_1    = "mono_stade_latest.geojson",
    fichier_2    = "mono_bivoltB_gen2_latest.geojson",
    nom_couche_1 = "Stade + bivoltisme A",
    nom_couche_2 = "Bivoltisme B (gen2)",
    champ_label_1 = "label",
    champ_label_2 = "label"
  )
)

ui <- fluidPage(
  titlePanel("Phenologie forestiere — France metropolitaine"),
  fluidRow(
    column(4,
      selectInput("modele", "Modele",
                  choices  = setNames(names(MODELES), sapply(MODELES, `[[`, "label")),
                  selected = "phenips")
    ),
    column(8,
      div(style = "margin-top: 25px; color: #666; font-size: 13px;",
          textOutput("statut", inline = TRUE))
    )
  ),
  leafletOutput("carte", height = "78vh")
)

server <- function(input, output, session) {

  lire_geojson_securise <- function(url) {
    tryCatch({
      fichier_tmp <- tempfile(fileext = ".geojson")
      request(url) |> req_perform(path = fichier_tmp)
      sf::st_read(fichier_tmp, quiet = TRUE)
    }, error = function(e) {
      message("Erreur lecture GeoJSON (", url, "): ", conditionMessage(e))
      NULL
    })
  }

  # Cache en memoire : evite de retelecharger si l'utilisateur revient sur
  # un modele deja consulte dans la meme session.
  cache_geojson <- reactiveValues()

  charger_couche <- function(cle_modele, suffixe, nom_fichier) {
    cle <- paste0(cle_modele, "_", suffixe)
    if (is.null(cache_geojson[[cle]])) {
      cache_geojson[[cle]] <- lire_geojson_securise(url_geojson(nom_fichier))
    }
    cache_geojson[[cle]]
  }

  couche1 <- reactive({
    m <- MODELES[[input$modele]]
    charger_couche(input$modele, "1", m$fichier_1)
  })

  couche2 <- reactive({
    m <- MODELES[[input$modele]]
    charger_couche(input$modele, "2", m$fichier_2)
  })

  output$statut <- renderText({
    m  <- MODELES[[input$modele]]
    n1 <- if (is.null(couche1())) 0 else nrow(couche1())
    n2 <- if (is.null(couche2())) 0 else nrow(couche2())
    glue("{m$label} — {n1} polygones ({m$nom_couche_1}) | {n2} polygones ({m$nom_couche_2})")
  })

  output$carte <- renderLeaflet({
    leaflet() |>
      addProviderTiles("CartoDB.Positron") |>
      setView(lng = 2.5, lat = 46.5, zoom = 6)
  })

  # --- Couche 1 : toujours un champ "couleur" pre-calcule cote export R ---
  observe({
    m <- MODELES[[input$modele]]
    g <- couche1()
    if (is.null(g) || nrow(g) == 0) return(invisible(NULL))

    couleurs <- if ("couleur" %in% names(g)) g$couleur else "#31688E"
    labels_1 <- if (m$champ_label_1 %in% names(g)) g[[m$champ_label_1]] else NULL

    leafletProxy("carte") |>
      clearGroup(m$nom_couche_1) |>
      addPolygons(
        data = g, fillColor = couleurs, color = "#333333", weight = 0.5,
        fillOpacity = 0.65, group = m$nom_couche_1, label = labels_1
      )
  })

  # --- Couche 2 : champ "couleur" pre-calcule, SAUF onset PHENIPS
  #     (date_iso continue -> degrade colorNumeric cote client) ---
  observe({
    m <- MODELES[[input$modele]]
    o <- couche2()
    if (is.null(o) || nrow(o) == 0) return(invisible(NULL))

    labels_2 <- if (m$champ_label_2 %in% names(o)) o[[m$champ_label_2]] else NULL

    if ("date_iso" %in% names(o)) {
      dates_num <- as.Date(o$date_iso)
      pal_onset <- colorNumeric("YlOrRd", domain = dates_num, na.color = "#cccccc")
      leafletProxy("carte") |>
        clearGroup(m$nom_couche_2) |>
        removeControl("legende_couche2") |>
        addPolygons(
          data = o, fillColor = ~pal_onset(as.Date(date_iso)), color = "#333333",
          weight = 0.5, fillOpacity = 0.75, group = m$nom_couche_2, label = labels_2
        ) |>
        addLegend(position = "bottomright", pal = pal_onset, values = dates_num,
                  title = "1er essaimage", labFormat = labelFormat(date = TRUE),
                  layerId = "legende_couche2", group = m$nom_couche_2)
    } else {
      couleurs <- if ("couleur" %in% names(o)) o$couleur else "#31688E"
      leafletProxy("carte") |>
        clearGroup(m$nom_couche_2) |>
        removeControl("legende_couche2") |>
        addPolygons(
          data = o, fillColor = couleurs, color = "#333333", weight = 0.5,
          fillOpacity = 0.65, group = m$nom_couche_2, label = labels_2
        )
    }
  })

  observe({
    m <- MODELES[[input$modele]]
    autres_groupes <- setdiff(
      unlist(lapply(MODELES, function(x) c(x$nom_couche_1, x$nom_couche_2))),
      c(m$nom_couche_1, m$nom_couche_2)
    )
    leafletProxy("carte") |>
      clearGroup(autres_groupes) |>
      addLayersControl(
        overlayGroups = c(m$nom_couche_1, m$nom_couche_2),
        options       = layersControlOptions(collapsed = FALSE)
      )
  })
}

shinyApp(ui, server)
