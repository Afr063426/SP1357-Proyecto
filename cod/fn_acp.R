# Se construye una funcion que retorna 
fn_acp <- function(df, scale = TRUE, name_acp, eje_x = "PC1",  eje_y = "PC2") {
    require(tidyverse) # To manipulate data

    require(ggforce) # To make the unit circle

    require(kableExtra) # To make beatiful tables

    # We select the name of columns of individuals
    name_observations <- colnames(df)[1]

    # We delete the first columns that have names
    names_df <- df[, 1]


    # We select the variables
    df <- df[, -1]


    # We get the number of observations
    n <- nrow(df)

    # We get the center of gravity
    g <- apply(df, 2, mean)

    # We estimate the standard deviation
    sd_n <- sqrt((n - 1) / n) * apply(df, 2, sd)

    # We assign a name to the acp
    df <- (df - matrix(rep(g, n), nrow = n, byrow = TRUE)) /
        matrix(rep(sd_n, n), nrow = n, byrow = TRUE)

    # We assign the name to the ACP
    assign(name_acp, prcomp(df), envir = .GlobalEnv)

    # We get the acp to use it
    acp <- get(name_acp)

    # We make the main plain
    df_pc <- as.data.frame(as.matrix(df) %*% as.matrix(acp$rotation))

    # We add the name of the observations
    df_pc <- cbind(names_df, df_pc)

    # We add the name of the first column
    colnames(df_pc)[1] <- name_observations

    # We select the names of the variables
    names_cols <- colnames(df)


    # We print the ggplot that is the main plain
    plot_plano_pc <- (ggplot(, aes(x = df_pc[,eje_x], y = df_pc[,eje_y], label = df_pc[, 1])) +
        geom_label(vjust = 1) +
        geom_point() +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 0) +
        theme_minimal() +
        # theme(text = element_text(size = 16)) +
        labs(x = eje_x, y = eje_y, caption = paste("Porcentaje de inercia:",
            sum((acp$sdev)[c(1, 2)]^(2)) / sum(acp$sdev^(2)),
            sep = " "
        )))
    
    # We estimate the correlation of the principal components
    df_correlations <- data.frame("Variable" = names_cols[1], (n - 1) /
        n * cor(df[, names_cols[1]], df_pc[, -1]))

    for (i in colnames(df)[-1]) {
        df_aux <- data.frame("Variable" = i, (n - 1) / n * cor(df[, i], df_pc[, -1]))
        df_correlations <- rbind(df_correlations, df_aux)
    }
    
    # We show the correlation circle of variables
    plot_c_correlaciones <- (ggplot(, aes(
        x = df_correlations[, eje_x], y = df_correlations[, eje_y],
        label = df_correlations$Variable
    )) +
        geom_label(vjust = 1) +
        geom_point() +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 0) +
        geom_circle(aes(x0 = 0, y0 = 0, r = 1)) +
        coord_fixed() +
        theme_minimal() +
        # theme(text = element_text(size = 16)) +
        labs(x = "PC1", y = "PC2", caption = paste("Porcentaje de inercia:",
            sum((acp$sdev)[c(1, 2)]^(2)) / sum(acp$sdev^(2)),
            sep = " "
        )))


    

    # We save the eigen values
    eigen_values <- data.frame(
        "Valor propio" = 1:ncol(df),
        "Valor" = acp$sdev^(2), "Porcentaje de inercia" = acp$sdev^(2) /
            sum(acp$sdev^(2)) * 100,
        "Porcentaje de inercia acumulado" = cumsum(acp$sdev^(2)) /
            sum(acp$sdev^(2)) * 100
    )

    # We show the table with values
    plot_elbow <- (ggplot(, aes(x = 1:length(acp$sdev), y = acp$sdev^(2))) +
        geom_line() +
        geom_point() +
        theme_minimal() +
        # theme(text = element_text(size = 16)) +
        labs(x = "Número", y = "Valor propio"))

    

    # We estimate the euclidean norm squared
    df_quality_representation <- apply(df, 1, function(x) sum((x)^(2)))

    # We estimate the quality of representation of each projection
    df_quality_representation <- df_pc[, -1]^2 /
        matrix(rep(df_quality_representation, ncol(df)), ncol = ncol(df))

    # We add the name of the columns
    df_quality_representation <- cbind(
        names_df,
        df_quality_representation
    )

    # We estimate the quality of the main plain
    quality_main_plain_id <- df_quality_representation[, 2] + df_quality_representation[, 3]

    df_quality_representation$Calidad.representacion.plano.principal <- quality_main_plain_id

    # We add the name of the column
    colnames(df_quality_representation)[1] <- name_observations

    #We estimate the communalities
    df_comunalities <- df_correlations %>%
        mutate("Comunalidad" = PC1^2 + PC2^2) %>%
        select(Variable, Comunalidad)
    
    #We return the list the list with the tables
    list(
        eigen_values = eigen_values,
        df_pc = df_pc,
        df_correlations = df_correlations,
        df_quality_representation = df_quality_representation,
        df_comunalities = df_comunalities,
        plot_plano_pc = plot_plano_pc,
        plot_elbow = plot_elbow, 
        plot_c_correlaciones = plot_c_correlaciones
    )
}

