packages <- c('reshape2','Cairo','tidyverse','reticulate','pracma','data.table',
              'grid',"plotly","shinyalert","shinydashboard","shinycssloaders",
              "rhandsontable","tableHTML")

invisible(lapply(packages, library, character.only = TRUE))

appName     <- "MoltenProt"
user        <- Sys.info()['user']

reticulate::use_python(paste0("/home/",user,"/myenv/bin/python"), required = TRUE)

# developer path
base_dir <- paste0("/home/",user,"/spc_shiny_servers/differentialScanningFluorimetryApps/appFiles/",appName,"/")

# path for the docker user
if (user == 'shiny') {
  base_dir <- paste0("/home/shiny/",appName,'/')
}

global_chunck_n     <- 16 # should be global_plot_columns * global_plot_rows
global_plot_columns <- 4
global_plot_rows    <- 4

myrenderer <- "function(instance, td, row, col, prop, value, cellProperties) {
                Handsontable.renderers.TextRenderer.apply(this, arguments);
                if (instance.params) {
                    hcols = instance.params.col_highlight
                    hcols = hcols instanceof Array ? hcols : [hcols]
                    hrows = instance.params.row_highlight
                    hrows = hrows instanceof Array ? hrows : [hrows]

                    for (i = 0; i < hcols.length; i++) {
                        if (hcols[i] == col && hrows[i] == row) {
                            td.style.background = instance.getDataAtCell(row, col);
                        }
                    }
                }
  }"

myrendererBoolean <- "function(instance, td, row, col, prop, value, cellProperties) {
            Handsontable.renderers.CheckboxRenderer.apply(this, arguments);}"

global_palette_9 <- c(
    "#E41A1C", "#377EB8", "#4DAF4A",
    "#984EA3", "#FF7F00", "#FFFF33",
    "#A65628", "#F781BF", "#999999"
)

global_palette_40 = c(
    "#8DD3C7", "#ADDFC1", "#CDEBBB", "#EDF8B6", "#F6F6B8", "#E4E2C3", "#D2CFCE",
    "#BFBBD9", "#CDABBF", "#DE9AA2", "#F08A84", "#EE857B", "#CB9297", "#A9A0B2",
    "#86AECE", "#9CB1B8", "#C0B299", "#E3B379", "#F7B762", "#E2C364", "#CDCE66",
    "#B8DA68", "#C1DA82", "#D6D5A5", "#EBD0C8", "#FACDE4", "#F0D1E1", "#E6D4DD",
    "#DCD7DA", "#D3C9D3", "#CBAFCC", "#C396C4", "#BC82BD", "#C0A0BF", "#C5BFC1",
    "#C9DDC3", "#D3EBB7", "#E2EB9F", "#F0EC87", "#FFED6F"
)

global_palette_96 <- c(
"#8DD3C7", "#9AD8C4", "#A7DDC2", "#B4E2C0", "#C1E7BD", "#CFECBB", "#DCF1B9", "#E9F6B6",
"#F6FBB4", "#FCFCB4", "#F4F4B9", "#EDECBD", "#E5E4C2", "#DEDCC6", "#D6D4CB", "#CFCCCF",
"#C7C4D4", "#C0BCD8", "#C3B5D1", "#CAAEC5", "#D1A7B9", "#D8A0AD", "#DF9AA1", "#E69395",
"#ED8C88", "#F4867C", "#F98073", "#EB867E", "#DD8B89", "#CE9194", "#C0979F", "#B29CAB",
"#A4A2B6", "#96A8C1", "#87ADCC", "#86B1CD", "#95B1BF", "#A3B1B2", "#B1B2A5", "#C0B298",
"#CEB28B", "#DDB37E", "#EBB371", "#FAB364", "#F5B762", "#EDBC63", "#E4C164", "#DCC665",
"#D3CB65", "#CBD066", "#C2D567", "#BADA68", "#B4DD6B", "#BCDB79", "#C5D988", "#CDD796",
"#D6D5A5", "#DED3B3", "#E7D1C1", "#EFCFD0", "#F8CDDE", "#F9CDE4", "#F5CFE2", "#F1D0E1",
"#EDD1E0", "#E9D3DE", "#E5D4DD", "#E1D6DB", "#DDD7DA", "#D9D8D9", "#D5CFD6", "#D2C5D2",
"#CFBBCF", "#CBB0CC", "#C8A6C9", "#C59CC5", "#C191C2", "#BE87BF", "#BC83BD", "#BE8FBE",
"#C09CBF", "#C2A8C0", "#C3B4C0", "#C5C1C1", "#C7CDC2", "#C9DAC3", "#CBE6C4", "#CFEBBE",
"#D5EBB4", "#DBEBAA", "#E1EBA0", "#E7EC96", "#EDEC8C", "#F3EC82", "#F9EC78", "#FFED6F"
)

global_palette_364 <- c(
  "#8DD3C7", "#90D4C6", "#93D5C5", "#97D7C5", "#9AD8C4", "#9ED9C3", "#A1DBC3",
  "#A5DCC2", "#A8DDC2", "#ACDFC1", "#AFE0C0", "#B3E1C0", "#B6E3BF", "#B9E4BF",
  "#BDE5BE", "#C0E7BD", "#C4E8BD", "#C7E9BC", "#CBEBBC", "#CEECBB", "#D2EDBA",
  "#D5EFBA", "#D9F0B9", "#DCF1B9", "#DFF3B8", "#E3F4B7", "#E6F5B7", "#EAF7B6",
  "#EDF8B6", "#F1F9B5", "#F4FBB4", "#F8FCB4", "#FBFDB3", "#FFFFB3", "#FDFCB4",
  "#FBFAB5", "#F9F8B6", "#F7F6B7", "#F5F4B8", "#F3F2BA", "#F1F0BB", "#EFEEBC",
  "#EDECBD", "#EBEABE", "#E9E8C0", "#E7E5C1", "#E5E3C2", "#E3E1C3", "#E1DFC4",
  "#DFDDC5", "#DDDBC7", "#DBD9C8", "#D9D7C9", "#D7D5CA", "#D5D3CB", "#D3D1CC",
  "#D1CECE", "#CFCCCF", "#CDCAD0", "#CBC8D1", "#C9C6D2", "#C7C4D4", "#C5C2D5",
  "#C3C0D6", "#C1BED7", "#BFBCD8", "#BEBADA", "#BFB8D6", "#C1B6D3", "#C3B4D0",
  "#C5B2CD", "#C7B1CA", "#C9AFC7", "#CAADC3", "#CCABC0", "#CEAABD", "#D0A8BA",
  "#D2A6B7", "#D4A4B4", "#D6A3B1", "#D7A1AD", "#D99FAA", "#DB9DA7", "#DD9CA4",
  "#DF9AA1", "#E1989E", "#E2969A", "#E49597", "#E69394", "#E89191", "#EA8F8E",
  "#EC8E8B", "#EE8C88", "#EF8A84", "#F18881", "#F3877E", "#F5857B", "#F78378",
  "#F98175", "#FB8072", "#F78174", "#F38277", "#EF847A", "#EC857D", "#E88780",
  "#E48883", "#E08A86", "#DD8B89", "#D98D8C", "#D58E8F", "#D19092", "#CE9195",
  "#CA9398", "#C6949B", "#C3969E", "#BF97A1", "#BB99A3", "#B79AA6", "#B49CA9",
  "#B09DAC", "#AC9FAF", "#A8A0B2", "#A5A2B5", "#A1A3B8", "#9DA5BB", "#9AA6BE",
  "#96A8C1", "#92A9C4", "#8EABC7", "#8BACCA", "#87AECD", "#83AFD0", "#80B1D3",
  "#83B1CF", "#87B1CC", "#8BB1C8", "#8FB1C5", "#92B1C1", "#96B1BE", "#9AB1BB",
  "#9EB1B7", "#A2B1B4", "#A5B1B0", "#A9B2AD", "#ADB2A9", "#B1B2A6", "#B5B2A3",
  "#B8B29F", "#BCB29C", "#C0B298", "#C4B295", "#C7B291", "#CBB28E", "#CFB28B",
  "#D3B387", "#D7B384", "#DAB380", "#DEB37D", "#E2B379", "#E6B376", "#EAB373",
  "#EDB36F", "#F1B36C", "#F5B368", "#F9B365", "#FDB462", "#FAB562", "#F8B662",
  "#F6B762", "#F4B962", "#F1BA63", "#EFBB63", "#EDBC63", "#EBBE63", "#E8BF63",
  "#E6C064", "#E4C264", "#E2C364", "#DFC464", "#DDC564", "#DBC765", "#D9C865",
  "#D6C965", "#D4CA65", "#D2CC66", "#D0CD66", "#CDCE66", "#CBD066", "#C9D166",
  "#C7D267", "#C4D367", "#C2D567", "#C0D667", "#BED767", "#BBD868", "#B9DA68",
  "#B7DB68", "#B5DC68", "#B3DE69", "#B5DD6C", "#B7DC70", "#B9DC74", "#BBDB78",
  "#BEDB7B", "#C0DA7F", "#C2DA83", "#C4D987", "#C6D98A", "#C9D88E", "#CBD892",
  "#CDD796", "#CFD799", "#D1D69D", "#D4D6A1", "#D6D5A5", "#D8D5A8", "#DAD4AC",
  "#DDD4B0", "#DFD3B4", "#E1D3B7", "#E3D2BB", "#E5D2BF", "#E8D1C3", "#EAD1C6",
  "#ECD0CA", "#EED0CE", "#F0CFD2", "#F3CFD5", "#F5CED9", "#F7CEDD", "#F9CDE1",
  "#FCCDE5", "#FACDE4", "#F9CDE4", "#F8CEE3", "#F7CEE3", "#F6CEE3", "#F5CFE2",
  "#F4CFE2", "#F3CFE2", "#F2D0E1", "#F1D0E1", "#F0D1E1", "#EFD1E0", "#EED1E0",
  "#EDD2DF", "#ECD2DF", "#EBD2DF", "#E9D3DE", "#E8D3DE", "#E7D3DE", "#E6D4DD",
  "#E5D4DD", "#E4D5DD", "#E3D5DC", "#E2D5DC", "#E1D6DB", "#E0D6DB", "#DFD6DB",
  "#DED7DA", "#DDD7DA", "#DCD7DA", "#DBD8D9", "#DAD8D9", "#D9D9D9", "#D8D6D8",
  "#D7D3D7", "#D6D0D6", "#D5CED5", "#D4CBD4", "#D3C8D3", "#D2C6D3", "#D1C3D2",
  "#D1C0D1", "#D0BED0", "#CFBBCF", "#CEB8CE", "#CDB5CD", "#CCB3CD", "#CBB0CC",
  "#CAADCB", "#CAABCA", "#C9A8C9", "#C8A5C8", "#C7A3C8", "#C6A0C7", "#C59DC6",
  "#C49AC5", "#C398C4", "#C395C3", "#C292C2", "#C190C2", "#C08DC1", "#BF8AC0",
  "#BE88BF", "#BD85BE", "#BC82BD", "#BC80BD", "#BC83BD", "#BC86BD", "#BD89BD",
  "#BD8CBD", "#BE90BE", "#BE93BE", "#BF96BE", "#BF99BE", "#C09DBF", "#C0A0BF",
  "#C1A3BF", "#C1A6BF", "#C2AAC0", "#C2ADC0", "#C3B0C0", "#C3B3C0", "#C4B7C1",
  "#C4BAC1", "#C5BDC1", "#C5C0C1", "#C6C4C2", "#C6C7C2", "#C7CAC2", "#C7CDC2",
  "#C8D1C3", "#C8D4C3", "#C9D7C3", "#C9DAC3", "#CADEC4", "#CAE1C4", "#CBE4C4",
  "#CBE7C4", "#CCEBC5", "#CDEBC2", "#CFEBBF", "#D0EBBD", "#D2EBBA", "#D3EBB7",
  "#D5EBB5", "#D6EBB2", "#D8EBB0", "#D9EBAD", "#DBEBAA", "#DDEBA8", "#DEEBA5",
  "#E0EBA3", "#E1EBA0", "#E3EB9D", "#E4EB9B", "#E6EC98", "#E7EC96", "#E9EC93",
  "#EAEC90", "#ECEC8E", "#EEEC8B", "#EFEC89", "#F1EC86", "#F2EC83", "#F4EC81",
  "#F5EC7E", "#F7EC7C", "#F8EC79", "#FAEC76", "#FBEC74", "#FDEC71", "#FFED6F"
)