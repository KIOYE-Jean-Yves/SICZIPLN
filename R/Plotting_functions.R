# fonction pour faire les graphiques
coef_plot_barre= function(matrice_coefficient, nom_axes=c("Columns of Y","Variables","Coefficients value","SICPLN"),grad_echelle=c(-5,5)){
  grad_max=grad_echelle[2]#max(matrice_coefficient)
  grad_min=grad_echelle[1]#min(matrice_coefficient)
  nom_axe_y=nom_axes[1]
  nom_axe_x=nom_axes[2]
  nom_gradient_couleur=nom_axes[3]
  titre=nom_axes[4]
  metlcoef_sic_genus=melt(matrice_coefficient)
  metlcoef_sic_genus$sparsity=(ifelse(metlcoef_sic_genus$value== 0, 'Zero', "Nonzero"))
  plot_genus2_sicpln=ggplot(data=metlcoef_sic_genus, aes(x=as.factor(Var2), y=as.factor(Var1),
                                                         pattern = sparsity, fill= value))+
    geom_tile_pattern(
      width = 0.9,          # valeur entre 0 et 1 (1 = pleine largeur)
      height = 0.9,         # valeur entre 0 et 1 (1 = pleine hauteur)pattern_color = NA,
      pattern_fill = "black",
      pattern_angle = 45,
      pattern_density = 0.015,
      pattern_spacing = 0.06,
      pattern_key_scale_factor = 1) +
    scale_pattern_manual(values = c(Zero = "circle", Nonzero = "none"),name="") +
    scale_fill_gradient2(low = "darkred", high = "darkgreen",mid="white", midpoint = 0, limit = c(grad_min,grad_max),name=nom_gradient_couleur,guide = guide_colourbar(
      barwidth = 0.3,   # largeur de la barre de légende (horizontal si légende verticale)
      barheight = 2.5     # hauteur totale de la légende
    ))+
    coord_equal() +
    labs(x = nom_axe_y,y = nom_axe_x,title = titre) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")))+theme(plot.title = element_text(size = 12, hjust = 0.5), axis.title = element_text(size = 12,colour = "black"), legend.text = element_text(size = 8,colour = "black"),legend.title = element_text(color = "black", size = 8),
                                                                              axis.text.x = element_text(size=5,colour = "black",angle = 45,hjust = 1),strip.text.x = element_text(size = 5, colour = "black"),axis.text.y = element_text(size = 8 ,colour = "black"))
  return(plot_genus2_sicpln)
}

coef_plot_barre_divise_species= function(matrice_coefficient, nom_axes=c("Phylum","Variables","Coefficients value","SICPLN"),grad_echelle=c(-5,5)){
  grad_max=grad_echelle[2]#max(matrice_coefficient)
  grad_min=grad_echelle[1]#min(matrice_coefficient)
  nom_axe_y=nom_axes[1]
  nom_axe_x=nom_axes[2]
  nom_gradient_couleur=nom_axes[3]
  titre=nom_axes[4]
  metlcoef_sic_genus=melt(matrice_coefficient)
  metlcoef_sic_genus$sparsity=(ifelse(metlcoef_sic_genus$value== 0, 'Zero', "Nonzero"))
  plot_genus2_sicpln=ggplot(data=metlcoef_sic_genus, aes(x=as.factor(Var2), y=as.factor(Var1),
                                                         pattern = sparsity, fill= value))+
    geom_tile_pattern(pattern_color = NA,
                      pattern_fill = "black",
                      pattern_angle = 45,
                      pattern_density = 0.15,
                      pattern_spacing = 0.015,
                      pattern_key_scale_factor = 1) +
    scale_pattern_manual(values = c(Zero = "circle", Nonzero = "none"),name="") +
    scale_fill_gradient2(low = "darkred", high = "darkgreen",mid="white", midpoint = 0, limit = c(grad_min,grad_max),name=nom_gradient_couleur)+
    coord_equal() +
    labs(x = nom_axe_y,y = nom_axe_x,title = titre) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")))+theme(plot.title = element_text(size = 12, hjust = 0.5), axis.title = element_text(size = 15,colour = "black"), legend.text = element_text(size = 15,colour = "black"),legend.title = element_text(color = "black", size = 15),
                                                                              axis.text.x = element_text(size=15,colour = "black",angle = 45,hjust = 1),strip.text.x = element_text(size = 10, colour = "black"),axis.text.y = element_text(size = 12 ,colour = "black"))
  return(plot_genus2_sicpln)
}
