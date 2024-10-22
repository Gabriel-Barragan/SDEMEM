Datos_naranja <- function(Plot=TRUE) {
  # Generar datos
  data("Orange")
  ggplot2::theme_set(theme_bw())
  
  # Estandarizar la edad
  # Orange$age <- Orange$age/max(Orange$age)

  summary(Orange)
  
  # Reordenar sujetos
  Orange$Arbol <- factor(Orange$Tree, 
                        level = 1:5)

  M <- Orange$Tree %>% unique(.) %>% length(.)
  No_obs <- Orange$age[Orange$Tree == 1] %>% length(.)
  
  datos <- array(0, dim=c(No_obs, 2, M))
  
  for (m in 1:M) {
    datos[,1,m] <- Orange$age[Orange$Tree == m]#/max(Orange$age[Orange$Tree == m])
    datos[,2,m] <- Orange$circumference[Orange$Tree == m]
  }
  
  # Diagrama de spaghetti - circunferencia
  if (Plot==TRUE) {
    p <- Orange %>%
      ggplot(., aes(x = age, y = circumference, color = Arbol)) +
      geom_line(aes(group = Tree)) + geom_point(size=4) +
      labs(x = "Edad (días)", y = "Circunferencia (mm)") +
      geom_line(linewidth=2) + theme(text = element_text(size = 20))  
  } else {
    p <- NULL
  }
  
  return(list(datos = datos, M = M, p = p))
}
