using PlotThemes, Colors, PlotUtils

let
mytheme_palette = [
    colorant"#d77255",
    colorant"#009afa",
    colorant"#707070",
    colorant"#21ab74",
    colorant"#ba3030",
    colorant"#9467bd"
]

_mytheme = PlotTheme(Dict([
    :background => :white,
    :framestyle => :box,
    :grid => false,
    :minorticks => 5,
    :linewidth => 1.4,
    :markerstrokewidth => 0,
    :fontfamily => "Computer Modern",
    :colorgradient => :magma,
    :guidefontsize => 12,
    :titlefontsize => 12,
    :tickfontsize => 8,
    :palette => mytheme_palette,
    :legend => :outertopright])
)

add_theme(:mytheme, _mytheme)
end
