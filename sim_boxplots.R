file21 <- readRDS("/Users/jingyihao/Desktop/harm_prediction_project/00_result/sim200_new/simulation200_1_250.RDS")
file22 <- readRDS("/Users/jingyihao/Desktop/harm_prediction_project/00_result/sim200_new/simulation200_251_500.RDS")
file23 <- readRDS("/Users/jingyihao/Desktop/harm_prediction_project/00_result/sim200_new/simulation200_501_750.RDS")
file24 <- readRDS("/Users/jingyihao/Desktop/harm_prediction_project/00_result/sim200_new/simulation200_751_1050.RDS")
file51 <- readRDS("/Users/jingyihao/Desktop/harm_prediction_project/00_result/sim500_new/simulation500_1_200.RDS")
file52 <- readRDS("/Users/jingyihao/Desktop/harm_prediction_project/00_result/sim500_new/simulation500_201_500.RDS")
file53 <- readRDS("/Users/jingyihao/Desktop/harm_prediction_project/00_result/sim500_new/simulation500_501_700.RDS")
file54 <- readRDS("/Users/jingyihao/Desktop/harm_prediction_project/00_result/sim500_new/simulation500_701_1050.RDS")

combined_data_200 <- rbind(file21, file22, file23, file24)
combined_data_500 <- rbind(file51, file52, file53, file54)


column_names <- c(paste("bias", 1:7, sep = ""), 
                  paste("se", 1:7, sep = ""), 
                  paste("mse", 1:7, sep = ""), 
                  paste("coverage", 1:7, sep = ""),
                  paste("convergence","_",sep=""),
                  paste("saddle","_",sep=""),
                  paste("PI_rate11","_",sep=""),
                  paste("PI_rate12","_",sep=""),
                  paste("PI_rate21","_",sep=""),
                  paste("PI_rate22","_",sep=""),
                  paste("PI_rate31","_",sep=""),
                  paste("PI_rate32","_",sep=""),
                  paste("PI_rate41","_",sep=""),
                  paste("PI_rate42","_",sep=""),
                  paste("PI_rate51","_",sep="")
                  )

# Assign column names to the data frame
colnames(combined_data_200) <- column_names
colnames(combined_data_500) <- column_names

# Step 1: Filter rows where the 29th column equals 1
filtered_data_200 <- combined_data_200[combined_data_200[,29] == 1, ]
filtered_data_500 <- combined_data_500[combined_data_500[,29] == 1, ]
# Step 2: Select the first 1000 rows from the filtered data
filtered_200 <- filtered_data_200[1:1000, ]
filtered_500 <- filtered_data_500[1:1000, ]

# Step 3: Calculate column means, removing NA values if present
final_200 <- colMeans(filtered_200)
final_500 <- colMeans(filtered_500)

colors <- c("lightgrey", "lightblue")
box_colors <- c(colors[1], colors[2], colors[1], colors[2], colors[1], 
                colors[2], colors[1], colors[2], colors[2])

short_labels <- c(expression(p[i1]), 
                  expression(p[i2]), 
                  expression(p[i1] + p[i2]), 
                  expression(p[i2]/p[i1]), 
                  expression(n[i2]))
label_positions <- c(1.5, 4.5, 7.5, 10.5, 13)


first_row_labels <- c(expression("Condition on " * X[i]),
                      expression("Condition on " * n[i1] * ", " * n[i2] * ", " * n[i2] * ", " * X[i]),
                      expression("Condition on " * X[i]),
                      expression("Condition on " * n[i1] * ", " * n[i2] * ", " * n[i2] * ", " * X[i]),
                      expression("Condition on " * X[i]),
                      expression("Condition on " * n[i1] * ", " * n[i2] * ", " * n[i2] * ", " * X[i]),
                      expression("Condition on " * X[i]),
                      expression("Condition on " * n[i1] * ", " * n[i2] * ", " * n[i2] * ", " * X[i]),
                      expression("Condition on " * n[i1] * ", " * X[i]))


# sample size = 200
group_data_200 <- filtered_200[, 31:39]

# Set positions for boxes so that pairs stay closer
positions <- c(1, 2, 4, 5, 7, 8, 10, 11, 13)

# Create the box plot with adjusted positions and custom colors
boxplot(group_data_200, 
        at = positions,  # Adjust positions for closer pairs
        main = "Boxplot of PICP with 200 sample size", 
        ylab = "Percentage", 
        col = box_colors,  # Use the custom colors for each box
        names = rep("", length(positions)),  # Use short labels for the x-axis
        las = 1,               # Keep labels horizontal
        cex.axis = 1,        # Shrink the size of axis labels
        cex.lab = 1,         # Shrink the size of y-axis label
        cex.main = 1,
        ylim = c(0.65,1))      
mtext(side = 1, at = label_positions, text = short_labels, line = 1, cex = 1)
# Add a horizontal line at y = 0.95
abline(h = 0.95, col = "red", lwd = 2)
axis(2, at = 0.95, labels = "0.95", las = 1,cex.axis = 1)

# Add a legend to explain the colors
legend("bottomleft", 
       legend <- c(expression("Condition on " * X[i]), 
                  expression("Condition on " * n[i1] * ", " * n[i2] * ", " * n[i2] * ", " * X[i])),
       fill = colors, cex = 1,  bty = "n")  # Add a legend with no border

# samlpe size = 500
group_data_500 <- filtered_500[, 31:39]

# Create the box plot with adjusted positions and custom colors
boxplot(group_data_500, 
        at = positions,  # Adjust positions for closer pairs
        main = "Boxplot of PICP with 500 sample size", 
        ylab = "Percentage", 
        col = box_colors,  # Use the custom colors for each box
        names = rep("", length(positions)),  # Use short labels for the x-axis
        las = 1,               # Keep labels horizontal
        cex.axis = 1,        # Shrink the size of axis labels
        cex.lab = 1,         # Shrink the size of y-axis label
        cex.main = 1,
        ylim = c(0.65,1))      
mtext(side = 1, at = label_positions, text = short_labels, line = 1, cex = 1)
# Add a horizontal line at y = 0.95
abline(h = 0.95, col = "red", lwd = 2)
axis(2, at = 0.95, labels = "0.95", las = 1,cex.axis = 1)

# Add a legend to explain the colors
legend("bottomleft", 
       legend <- c(expression("Condition on " * X[i]), 
                  expression("Condition on " * n[i1] * ", " * n[i2] * ", " * n[i2] * ", " * X[i])),
       fill = colors, cex = 1,  bty = "n")  # Add a legend with no border

png("combined_plots_202.png", width = 800, height = 300)  # Adjust width and height as needed
# Define the layout matrix and specify relative widths and heights
layout(matrix(c(1, 2), nrow = 1, byrow = TRUE), 
       widths = c(4, 4),   # Set relative widths for each column
       heights = c(3, 3))  # Set relative heights for each row
par(mar = c(3, 3, 3, 3))  # Adjust margins to allow space for titles

# Plot for sample size 200 with adjusted title
boxplot(group_data_200, 
        at = positions, 
        ylab = "Percentage", 
        col = box_colors, 
        names = rep("", length(positions)), 
        las = 1, 
        cex.axis = 1, 
        cex.lab = 1, 
        cex.main = 1, 
        ylim = c(0.65, 1)) 
mtext(side = 3, line = 1, adj =0, "(a) Boxplot of PICP with 200 sample size")
mtext(side = 1, at = label_positions, text = short_labels, line = 1, cex = 1)
abline(h = 0.95, col = "red", lwd = 2)
axis(2, at = 0.95, labels = "0.95", las = 1, cex.axis = 1)
legend("bottomleft", 
       legend <- c(expression("Condition on " * X[i]), 
                  expression("Condition on " * n[i1] * ", " * n[i2] * ", " * n[i2] * ", " * X[i])),
       fill = colors, cex = 1, bty = "n")

# Plot for sample size 500 with adjusted title
boxplot(group_data_500, 
        at = positions, 
        ylab = "Percentage", 
        col = box_colors, 
        names = rep("", length(positions)), 
        las = 1, 
        cex.axis = 1, 
        cex.lab = 1, 
        cex.main = 1, 
        ylim = c(0.65, 1))      
mtext(side = 3, line = 1, adj =0, "(b) Boxplot of PICP with 500 sample size")
mtext(side = 1, at = label_positions, text = short_labels, line = 1, cex = 1)
abline(h = 0.95, col = "red", lwd = 2)
axis(2, at = 0.95, labels = "0.95", las = 1, cex.axis = 1)
legend("bottomleft", 
       legend <- c(expression("Condition on " * X[i]), 
                  expression("Condition on " * n[i1] * ", " * n[i2] * ", " * n[i2] * ", " * X[i])),
       fill = colors, cex = 1, bty = "n")

# Reset the layout after plotting
par(mfrow = c(1,1))

dev.off()
