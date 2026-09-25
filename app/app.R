# EDGE2 x VGP dashboard
# Runs as a normal Shiny app (shiny::runApp("app")) or in the browser via
# shinylive (shinylive::export("app", "_site")). House style: Arial + Zissou1.
library(shiny)
library(bslib)
library(plotly)
library(DT)

# ---------------------------------------------------------------- data ----
sp   <- read.csv("data/species.csv", stringsAsFactors = FALSE, na.strings = c("NA"))
summ <- read.csv("data/summary.csv", stringsAsFactors = FALSE)
enr  <- read.csv("data/enrichment.csv", stringsAsFactors = FALSE)
plc  <- read.csv("data/vgp_placement.csv", stringsAsFactors = FALSE, na.strings = c("NA"))
sp$has_vgp_genome  <- as.logical(sp$has_vgp_genome)
sp$is_EDGE_species <- as.integer(sp$is_EDGE_species)

CLADES <- c(Mammals = "mammals", Birds = "birds", Squamates = "squamates",
            Amphibians = "amphibians", Chondrichthyans = "chondrichthyans")
CLADE_LAB <- setNames(names(CLADES), CLADES)
GROUP_LEVEL <- c(mammals = "order", birds = "order", squamates = "family",
                 amphibians = "order", chondrichthyans = "order")

ZISSOU <- c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")
RL_ORDER <- c("LC", "NT", "VU", "EN", "CR", "EW", "EX", "DD", "NM")
RL_COL <- c(LC = "#3B9AB2", NT = "#78B7C5", VU = "#EBCC2A", EN = "#E96500", CR = "#F21A00",
            EW = "#4A4A4A", EX = "#1A1A1A", DD = "#B8B8B8", NM = "#E6E6E6")
RL_NAME <- c(LC = "Least Concern", NT = "Near Threatened", VU = "Vulnerable",
             EN = "Endangered", CR = "Critically Endangered", EW = "Extinct in the Wild",
             EX = "Extinct", DD = "Data Deficient", NM = "No Red List match")
CLADE_COL <- c(mammals = "#3B9AB2", birds = "#EBCC2A", squamates = "#E96500",
               amphibians = "#F21A00", chondrichthyans = "#78B7C5")
FONT <- list(family = "Arial, Helvetica, sans-serif", size = 12, color = "#2B2B2B")

sp$RLcat <- factor(sp$RLcat, levels = RL_ORDER)
sp$label <- ifelse(nzchar(sp$common_name), paste0(sp$species, " (", sp$common_name, ")"), sp$species)
sp$genome <- ifelse(sp$has_vgp_genome, "VGP genome", "no genome")

hover_txt <- function(d) {
  paste0("<i>", d$species, "</i>",
         ifelse(nzchar(d$common_name), paste0("<br>", d$common_name), ""),
         "<br>", d$order, " / ", d$family,
         "<br>Red List: ", RL_NAME[as.character(d$RLcat)],
         "<br>EDGE2 rank #", d$EDGErank, " · ", sprintf("%.2f", d$EDGEmed), " Myr",
         "<br>ED2 ", sprintf("%.1f", d$EDmed), " Myr · pext ", sprintf("%.3f", d$pextmed),
         ifelse(d$is_EDGE_species == 1, "<br><b>EDGE species</b>", ""),
         ifelse(d$has_vgp_genome, paste0("<br>VGP genome: ", d$vgp_accession), "<br>no VGP genome"))
}

plotly_style <- function(p, ...) {
  p %>% layout(font = FONT, paper_bgcolor = "white", plot_bgcolor = "white",
               legend = list(orientation = "v", font = list(size = 11)), ...) %>%
    config(displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "select2d"),
           toImageButtonOptions = list(format = "svg", filename = "edge2_plot"))
}

dt_opts <- list(pageLength = 15, scrollX = TRUE, dom = "Bfrtip",
                buttons = list("copy", list(extend = "csv", filename = "edge2_selection")))

theme <- bs_theme(version = 5, primary = "#3B9AB2", secondary = "#78B7C5",
                  base_font = "Arial, Helvetica, sans-serif",
                  heading_font = "Arial, Helvetica, sans-serif", "navbar-bg" = "#3B9AB2")

vbox <- function(title, value, sub = NULL, col = "#3B9AB2") {
  value_box(title = title, value = value, p(sub), theme = value_box_theme(bg = col, fg = "white"))
}

# ------------------------------------------------------------------ UI ----
ui <- page_navbar(
  title = "EDGE2 × VGP", theme = theme, fillable = FALSE,
  nav_panel("Overview",
    layout_columns(col_widths = c(3, 3, 3, 3),
      vbox("Species scored", format(nrow(sp), big.mark = ","), "5 VertLife trees × 1,000 posterior trees"),
      vbox("EDGE species", format(sum(sp$is_EDGE_species), big.mark = ","), "threatened, ED2 ≥ median in ≥ 50% of trees", "#78B7C5"),
      vbox("EDGE species with a VGP genome", sum(sp$is_EDGE_species == 1 & sp$has_vgp_genome),
           sprintf("%.1f%% of all EDGE species", 100 * mean(sp$has_vgp_genome[sp$is_EDGE_species == 1])), "#E1AF00"),
      vbox("VGP genomes placed", format(sum(sp$has_vgp_genome), big.mark = ","), "on tree tips across the five groups", "#F21A00")),
    layout_columns(col_widths = c(6, 6),
      card(card_header("VGP coverage of EDGE species"), plotlyOutput("ov_cov", height = 340)),
      card(card_header("Coverage rises with distinctness and EDGE rank"), plotlyOutput("ov_enr", height = 340))),
    card(card_header("Summary by group"), DTOutput("ov_table"))
  ),
  nav_panel("Explore a group",
    layout_sidebar(
      sidebar = sidebar(width = 300,
        selectInput("clade", "Group", CLADES),
        checkboxGroupInput("rl", "Red List category", choices = setNames(RL_ORDER, paste(RL_ORDER, "–", RL_NAME[RL_ORDER])),
                           selected = RL_ORDER),
        selectizeInput("grp", "Order / family", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "all")),
        radioButtons("edge_only", "Species", c("All species" = "all", "EDGE species only" = "edge"), inline = FALSE),
        radioButtons("genome", "VGP genome", c("Any" = "any", "Has genome" = "yes", "No genome" = "no")),
        textInput("search", "Search name", placeholder = "e.g. pangolin, Manis"),
        helpText("Filters apply to every plot and the table on this tab.")),
      layout_columns(col_widths = c(3, 3, 3, 3),
        uiOutput("vb_n"), uiOutput("vb_edge"), uiOutput("vb_gen"), uiOutput("vb_edgegen")),
      navset_card_tab(
        nav_panel("Distinctness vs extinction risk", plotlyOutput("p_scatter", height = 520),
                  helpText("Each point is a species (median across 1,000 trees). Horizontal bands reflect the Red List category behind each extinction probability.")),
        nav_panel("Top species",
                  layout_columns(col_widths = c(4, 4, 4),
                    radioButtons("top_by", NULL, c("EDGE2 score" = "EDGEmed", "ED2 (distinctness)" = "EDmed"), inline = TRUE),
                    sliderInput("top_n", "Number of species", 10, 100, 40, step = 5),
                    div()),
                  plotlyOutput("p_top", height = 760)),
        nav_panel("EDGE2 rank curve", plotlyOutput("p_rank", height = 460)),
        nav_panel("By order / family", plotlyOutput("p_group", height = 560))
      ),
      card(card_header("Species table (filtered)"), DTOutput("t_species"))
    )
  ),
  nav_panel("Sequenced EDGE species",
    layout_sidebar(
      sidebar = sidebar(width = 260, checkboxGroupInput("seq_clades", "Groups", CLADES, selected = CLADES)),
      card(card_header("Every EDGE species that already has a VGP reference genome"),
           plotlyOutput("p_seq", height = 900)),
      card(DTOutput("t_seq")))
  ),
  nav_panel("Sequencing priorities",
    layout_sidebar(
      sidebar = sidebar(width = 280,
        selectInput("pr_clade", "Group", CLADES),
        selectizeInput("pr_grp", "Order / family", choices = NULL, multiple = TRUE, options = list(placeholder = "all")),
        checkboxGroupInput("pr_rl", "Red List category", c("VU", "EN", "CR"), selected = c("VU", "EN", "CR"), inline = TRUE),
        sliderInput("pr_n", "Show top", 10, 100, 30, step = 5)),
      card(card_header("Highest-ranked EDGE species without a VGP genome"), plotlyOutput("p_pr", height = 720)),
      card(DTOutput("t_pr")))
  ),
  nav_panel("VGP placement",
    card(card_header("How each VGP species was placed on a tree tip"),
         p("Every name in the VGP Ordinal List (Phase 1+ and Families sheets), with its placement tier, the tree tip it was matched to, and the evidence used. ",
           "Genomes in lineages without a VertLife tree (ray-finned fishes, turtles, crocodilians and others) are listed as out of scope."),
         DTOutput("t_plc")))
  ,
  nav_panel("About",
    card(
      h4("EDGE2 × VGP"),
      p("EDGE2 (Evolutionarily Distinct and Globally Endangered; Gumbs et al. 2023) rankings for five VertLife phylogenies, cross-referenced with the Vertebrate Genomes Project Ordinal List."),
      tags$ul(
        tags$li("Trees: Upham et al. 2019 (mammals), Jetz et al. 2012 (birds), Tonini et al. 2016 (squamates), Jetz & Pyron 2018 (amphibians), Stein et al. 2018 (chondrichthyans); 1,000 posterior trees each."),
        tags$li("Endangerment: IUCN Red List v2026-1. NM = no Red List match (a project code, not an IUCN category); treated like DD."),
        tags$li("An EDGE species is threatened (VU/EN/CR) and has ED2 at or above the tree median in at least 50% of trees."),
        tags$li("A species has a VGP genome if its VGP entry has a GCA_/GCF_ assembly accession.")),
      p("Full methods, data tables and figures: ", a("github.com/tanyalama/project_EDGE2", href = "https://github.com/tanyalama/project_EDGE2", target = "_blank"))))
)

# -------------------------------------------------------------- server ----
server <- function(input, output, session) {

  # ---- overview
  output$ov_cov <- renderPlotly({
    s <- summ[match(rev(CLADES), summ$clade), ]
    lab <- CLADE_LAB[s$clade]
    plot_ly(y = factor(lab, levels = lab)) %>%
      add_bars(x = ~s$n_EDGE_with_genome, name = "with VGP genome", orientation = "h", marker = list(color = "#3B9AB2"),
               hovertemplate = "%{x} with genome<extra></extra>") %>%
      add_bars(x = ~(s$n_EDGE - s$n_EDGE_with_genome), name = "without genome", orientation = "h", marker = list(color = "#D9D9D9"),
               text = sprintf("%d / %s (%.1f%%)", s$n_EDGE_with_genome, format(s$n_EDGE, big.mark = ","), s$pct_EDGE_with_genome),
               textposition = "outside", hovertemplate = "%{x} without genome<extra></extra>") %>%
      layout(barmode = "stack", xaxis = list(title = "EDGE species (n)"), yaxis = list(title = "")) %>% plotly_style()
  })
  output$ov_enr <- renderPlotly({
    e <- enr[match(rev(CLADES), enr$clade), ]; lab <- factor(CLADE_LAB[e$clade], levels = CLADE_LAB[e$clade])
    ser <- list(c("pct_species", "All species", "#B8B8B8", "circle"),
                c("pct_ED2", "Share of clade ED2 captured", "#78B7C5", "circle"),
                c("pct_EDGE_species", "EDGE species", "#FFFFFF", "square"),
                c("pct_top100_ED", "Top 100 species by ED2", "#3B9AB2", "diamond"),
                c("pct_top25_EDGE", "Top 25 EDGE species", "#2B2B2B", "triangle-up"))
    p <- plot_ly()
    for (s in ser) p <- p %>% add_markers(x = e[[s[1]]], y = lab, name = s[2],
        marker = list(color = s[3], symbol = s[4], size = 11, line = list(color = "#2B2B2B", width = 1)),
        hovertemplate = paste0(s[2], ": %{x:.1f}%<extra>%{y}</extra>"))
    p %>% layout(xaxis = list(title = "With a VGP reference genome (%)", rangemode = "tozero"), yaxis = list(title = "")) %>% plotly_style()
  })
  output$ov_table <- renderDT({
    s <- summ[match(CLADES, summ$clade), ]
    out <- data.frame(Group = CLADE_LAB[s$clade], Species = s$n_species, Threatened = s$n_threatened,
                      `EDGE species` = s$n_EDGE, `EDGE with VGP genome` = s$n_EDGE_with_genome,
                      `% EDGE with genome` = round(s$pct_EDGE_with_genome, 1),
                      `Top-50 EDGE with genome` = s$top50_EDGE_with_genome,
                      `VGP genomes placed` = s$n_vgp_genome,
                      `Expected PD loss (%)` = round(s$pct_PD_at_risk, 1),
                      `Expected PD loss (Gy)` = round(s$ePDloss_Gy_median, 1),
                      `No Red List match` = s$n_no_RedList_match, check.names = FALSE)
    datatable(out, rownames = FALSE, options = list(dom = "t", ordering = FALSE))
  })

  # ---- explore
  observeEvent(input$clade, {
    lv <- GROUP_LEVEL[[input$clade]]
    ch <- sort(unique(sp[sp$clade == input$clade, lv]))
    updateSelectizeInput(session, "grp", label = if (lv == "family") "Family" else "Order", choices = ch, selected = character(0), server = TRUE)
  })
  filt <- reactive({
    d <- sp[sp$clade == input$clade, ]
    d <- d[as.character(d$RLcat) %in% input$rl, ]
    if (length(input$grp)) d <- d[d[[GROUP_LEVEL[[input$clade]]]] %in% input$grp, ]
    if (input$edge_only == "edge") d <- d[d$is_EDGE_species == 1, ]
    if (input$genome == "yes") d <- d[d$has_vgp_genome, ]
    if (input$genome == "no")  d <- d[!d$has_vgp_genome, ]
    q <- trimws(input$search)
    if (nzchar(q)) d <- d[grepl(q, d$species, ignore.case = TRUE) | grepl(q, d$common_name, ignore.case = TRUE) |
                           grepl(q, d$order, ignore.case = TRUE) | grepl(q, d$family, ignore.case = TRUE), ]
    d
  })
  small_vb <- function(title, val, col) value_box(title = title, value = val, theme = value_box_theme(bg = col, fg = "white"), height = "110px")
  output$vb_n <- renderUI(small_vb("Species shown", format(nrow(filt()), big.mark = ","), "#3B9AB2"))
  output$vb_edge <- renderUI(small_vb("EDGE species", format(sum(filt()$is_EDGE_species), big.mark = ","), "#78B7C5"))
  output$vb_gen <- renderUI(small_vb("With VGP genome", sum(filt()$has_vgp_genome), "#E1AF00"))
  output$vb_edgegen <- renderUI(small_vb("EDGE with genome", sum(filt()$has_vgp_genome & filt()$is_EDGE_species == 1), "#F21A00"))

  output$p_scatter <- renderPlotly({
    d <- filt(); validate(need(nrow(d) > 0, "No species match the current filters."))
    set.seed(1); d$x <- d$pextmed * exp(runif(nrow(d), -0.08, 0.08))
    p <- plot_ly()
    for (cat in intersect(c("NM", "DD", "LC", "NT", "VU", "EN", "CR", "EW", "EX"), as.character(unique(d$RLcat)))) {
      s <- d[as.character(d$RLcat) == cat, ]
      p <- p %>% add_trace(type = "scattergl", mode = "markers", x = s$x, y = s$EDmed, name = RL_NAME[cat],
                           text = hover_txt(s), hoverinfo = "text",
                           marker = list(color = RL_COL[cat], size = ifelse(s$has_vgp_genome, 9, 5),
                                         line = list(color = ifelse(s$has_vgp_genome, "#2B2B2B", RL_COL[cat]), width = ifelse(s$has_vgp_genome, 1.2, 0)),
                                         opacity = 0.85))
    }
    p %>% layout(xaxis = list(title = "GE2: extinction probability (pext, median)", type = "log"),
                 yaxis = list(title = "ED2 (Myr, median)", type = "log"),
                 annotations = list(list(text = "Larger outlined points have a VGP genome", x = 0, y = 1.06, xref = "paper", yref = "paper",
                                         showarrow = FALSE, xanchor = "left", font = list(size = 11)))) %>% plotly_style()
  })

  output$p_top <- renderPlotly({
    d <- filt(); validate(need(nrow(d) > 0, "No species match the current filters."))
    v <- input$top_by; iq <- if (v == "EDGEmed") "EDGEiqr" else "EDiqr"
    d <- head(d[order(-d[[v]]), ], input$top_n); d <- d[rev(seq_len(nrow(d))), ]
    yl <- factor(d$label, levels = d$label)
    plot_ly() %>%
      add_segments(x = pmax(d[[v]] - d[[iq]] / 2, 0), xend = d[[v]] + d[[iq]] / 2, y = yl, yend = yl,
                   line = list(color = "#8C8C8C", width = 2), showlegend = FALSE, hoverinfo = "skip") %>%
      add_markers(x = d[[v]], y = yl, text = hover_txt(d), hoverinfo = "text", showlegend = FALSE,
                  marker = list(color = RL_COL[as.character(d$RLcat)], size = 11, line = list(color = "#2B2B2B", width = 1))) %>%
      add_markers(x = rep(max(d[[v]] + d[[iq]] / 2) * 1.06, nrow(d)), y = yl, name = "VGP genome", hoverinfo = "text",
                  text = ifelse(d$has_vgp_genome, paste("VGP genome:", d$vgp_accession), "no VGP genome"),
                  marker = list(symbol = "square", size = 10, color = ifelse(d$has_vgp_genome, "#3B9AB2", "white"),
                                line = list(color = "#3B9AB2", width = 1.5)), showlegend = FALSE) %>%
      layout(xaxis = list(title = if (v == "EDGEmed") "EDGE2 score (Myr; median ± ½ IQR)" else "ED2 (Myr; median ± ½ IQR)", rangemode = "tozero"),
             yaxis = list(title = "", tickfont = list(size = 10)),
             annotations = list(list(text = "■ = VGP genome", x = 1, y = 1.02, xref = "paper", yref = "paper", showarrow = FALSE,
                                     xanchor = "right", font = list(size = 11, color = "#3B9AB2")))) %>% plotly_style()
  })

  output$p_rank <- renderPlotly({
    d <- sp[sp$clade == input$clade, ]; d <- d[order(d$EDGErank), ]
    h <- filt()
    plot_ly() %>%
      add_ribbons(x = d$EDGErank, ymin = pmax(d$EDGEmed - d$EDGEiqr / 2, 0), ymax = d$EDGEmed + d$EDGEiqr / 2,
                  fillcolor = "rgba(120,183,197,0.45)", line = list(width = 0), name = "± ½ IQR across trees", hoverinfo = "skip") %>%
      add_lines(x = d$EDGErank, y = d$EDGEmed, line = list(color = "#3B9AB2", width = 2), name = "Median EDGE2", hoverinfo = "skip") %>%
      add_trace(type = "scattergl", mode = "markers", x = h$EDGErank, y = h$EDGEmed, name = "Filtered species",
                marker = list(color = RL_COL[as.character(h$RLcat)], size = 6, line = list(color = "#2B2B2B", width = 0.5)),
                text = hover_txt(h), hoverinfo = "text") %>%
      layout(xaxis = list(title = "EDGE2 rank (log scale)", type = "log"), yaxis = list(title = "EDGE2 score (Myr, median)", rangemode = "tozero")) %>%
      plotly_style()
  })

  output$p_group <- renderPlotly({
    d <- filt(); validate(need(nrow(d) > 0, "No species match the current filters."))
    lv <- GROUP_LEVEL[[input$clade]]; E <- d[d$is_EDGE_species == 1, ]
    validate(need(nrow(E) > 0, "No EDGE species in the current selection."))
    tab <- aggregate(cbind(n = 1, g = E$has_vgp_genome) ~ E[[lv]], FUN = sum)
    names(tab)[1] <- "grp"; tab <- head(tab[order(-tab$n), ], 25); tab <- tab[rev(seq_len(nrow(tab))), ]
    yl <- factor(tab$grp, levels = tab$grp)
    plot_ly(y = yl) %>%
      add_bars(x = tab$g, name = "with VGP genome", orientation = "h", marker = list(color = "#3B9AB2")) %>%
      add_bars(x = tab$n - tab$g, name = "without genome", orientation = "h", marker = list(color = "#D9D9D9"),
               text = paste0(tab$g, "/", tab$n), textposition = "outside") %>%
      layout(barmode = "stack", xaxis = list(title = "EDGE species (n)"), yaxis = list(title = ""),
             title = list(text = paste0("EDGE species by ", lv, " (top 25 in selection)"), font = list(size = 13), x = 0)) %>%
      plotly_style()
  })

  tbl_cols <- c("EDGErank", "species", "common_name", "order", "family", "RLcat", "EDGEmed", "EDmed", "pextmed",
                "isEDGEsp_frac", "is_EDGE_species", "has_vgp_genome", "vgp_accession", "vgp_placement_tier")
  tbl_names <- c("EDGE2 rank", "Species", "Common name", "Order", "Family", "Red List", "EDGE2 (Myr)", "ED2 (Myr)", "pext",
                 "Frac. trees EDGE", "EDGE species", "VGP genome", "VGP accession", "VGP placement")
  output$t_species <- renderDT({
    d <- filt()[, tbl_cols]; names(d) <- tbl_names
    datatable(d, rownames = FALSE, extensions = "Buttons", filter = "top", options = dt_opts) %>%
      formatStyle("Species", fontStyle = "italic") %>%
      formatStyle("Red List", backgroundColor = styleEqual(RL_ORDER, RL_COL[RL_ORDER]),
                  color = styleEqual(c("EX", "EW"), c("white", "white")))
  })

  # ---- sequenced EDGE species
  seqd <- reactive({
    d <- sp[sp$is_EDGE_species == 1 & sp$has_vgp_genome & sp$clade %in% input$seq_clades, ]
    d[order(match(d$clade, CLADES), d$EDGErank), ]
  })
  output$p_seq <- renderPlotly({
    d <- seqd(); validate(need(nrow(d) > 0, "Select at least one group."))
    d <- d[rev(seq_len(nrow(d))), ]
    yl <- factor(paste0(d$label, "  [", CLADE_LAB[d$clade], " #", d$EDGErank, "]"),
                 levels = paste0(d$label, "  [", CLADE_LAB[d$clade], " #", d$EDGErank, "]"))
    plot_ly() %>%
      add_segments(x = 0, xend = d$EDGEmed, y = yl, yend = yl, line = list(color = "#8C8C8C", width = 1.5), showlegend = FALSE, hoverinfo = "skip") %>%
      add_markers(x = d$EDGEmed, y = yl, text = hover_txt(d), hoverinfo = "text", showlegend = FALSE,
                  marker = list(color = RL_COL[as.character(d$RLcat)], size = 11, line = list(color = "#2B2B2B", width = 1))) %>%
      layout(xaxis = list(title = "EDGE2 score (Myr, median)", type = "log"), yaxis = list(title = "", tickfont = list(size = 10)),
             height = max(400, 18 * nrow(d) + 80)) %>% plotly_style()
  })
  output$t_seq <- renderDT({
    d <- seqd(); d$Group <- CLADE_LAB[d$clade]
    d <- d[, c("Group", "EDGErank", "species", "common_name", "order", "RLcat", "EDGEmed", "EDmed", "vgp_name", "vgp_accession", "vgp_placement_tier")]
    names(d) <- c("Group", "EDGE2 rank", "Species", "Common name", "Order", "Red List", "EDGE2 (Myr)", "ED2 (Myr)", "VGP name", "VGP accession", "VGP placement")
    datatable(d, rownames = FALSE, extensions = "Buttons", options = dt_opts) %>% formatStyle("Species", fontStyle = "italic")
  })

  # ---- priorities
  observeEvent(input$pr_clade, {
    lv <- GROUP_LEVEL[[input$pr_clade]]
    E <- sp[sp$clade == input$pr_clade & sp$is_EDGE_species == 1 & !sp$has_vgp_genome, ]
    updateSelectizeInput(session, "pr_grp", label = if (lv == "family") "Family" else "Order",
                         choices = sort(unique(E[[lv]])), selected = character(0), server = TRUE)
  })
  prd <- reactive({
    d <- sp[sp$clade == input$pr_clade & sp$is_EDGE_species == 1 & !sp$has_vgp_genome & as.character(sp$RLcat) %in% input$pr_rl, ]
    if (length(input$pr_grp)) d <- d[d[[GROUP_LEVEL[[input$pr_clade]]]] %in% input$pr_grp, ]
    d[order(d$EDGErank), ]
  })
  output$p_pr <- renderPlotly({
    d <- head(prd(), input$pr_n); validate(need(nrow(d) > 0, "No species match the current filters."))
    d <- d[rev(seq_len(nrow(d))), ]; yl <- factor(paste0(d$label, "  #", d$EDGErank), levels = paste0(d$label, "  #", d$EDGErank))
    plot_ly() %>%
      add_segments(x = pmax(d$EDGEmed - d$EDGEiqr / 2, 0), xend = d$EDGEmed + d$EDGEiqr / 2, y = yl, yend = yl,
                   line = list(color = "#8C8C8C", width = 2), showlegend = FALSE, hoverinfo = "skip") %>%
      add_markers(x = d$EDGEmed, y = yl, text = hover_txt(d), hoverinfo = "text", showlegend = FALSE,
                  marker = list(color = RL_COL[as.character(d$RLcat)], size = 11, line = list(color = "#2B2B2B", width = 1))) %>%
      layout(xaxis = list(title = "EDGE2 score (Myr; median ± ½ IQR)", rangemode = "tozero"), yaxis = list(title = "", tickfont = list(size = 10))) %>%
      plotly_style()
  })
  output$t_pr <- renderDT({
    d <- prd()[, c("EDGErank", "species", "common_name", "order", "family", "RLcat", "EDGEmed", "EDmed", "isEDGEsp_frac")]
    names(d) <- c("EDGE2 rank", "Species", "Common name", "Order", "Family", "Red List", "EDGE2 (Myr)", "ED2 (Myr)", "Frac. trees EDGE")
    datatable(d, rownames = FALSE, extensions = "Buttons", options = dt_opts) %>% formatStyle("Species", fontStyle = "italic")
  })

  # ---- placement
  output$t_plc <- renderDT({
    d <- plc[, c("vgp_name", "vgp_english_name", "vgp_lineage", "vgp_sheet", "vgp_accession", "clade", "placement_tier", "tree_tip", "credit_genome_to_tip", "evidence")]
    d$clade <- ifelse(is.na(d$clade), "", CLADE_LAB[d$clade])
    names(d) <- c("VGP name", "English name", "VGP lineage", "Sheet", "Accession", "Group", "Placement tier", "Tree tip", "Genome credited", "Evidence")
    datatable(d, rownames = FALSE, extensions = "Buttons", filter = "top", options = dt_opts)
  })
}

shinyApp(ui, server)
