import edu.princeton.cs.algs4.Picture;

public class SeamCarver {
    private static final double BORDER_ENERGY = 1000.00;
    private Picture picture;
    private int width;
    private int height;
        
    /** 
     * create a seam carver object based on the given picture
     * @param picture
     */

   public SeamCarver(Picture picture) {
       if (null == picture) {
           throw new IllegalArgumentException("Picture object is null");
       }
       
       this.picture = picture;
       width = this.picture.width();
       height = this.picture.height();
   }
   
    private double min(double a, double b, double c) {
        double minimum = a;
        if (minimum > b) {
            minimum = b;
        }
        
        if (minimum > c) {
            minimum = c;
        }
        
        return minimum;
    }
    
    private int minimumHIndex(double[][] cumEnergy, int minCol, int row) {
        
        int minIndex = minCol - 1;
        if (minCol - 1 < 0) {
            return cumEnergy[row][minCol] < cumEnergy[row][minCol + 1]
                    ? minCol : minCol + 1;
        }
        
        if (minCol + 1 > width - 1) {
            return cumEnergy[row][minCol - 1] < cumEnergy[row][minCol] 
                    ? minCol - 1 : minCol;
        }
        
        double minVal = cumEnergy[row][minCol - 1];
        if (minVal > cumEnergy[row][minCol]) {
            minVal = cumEnergy[row][minCol];
            minIndex = minCol;
        }
        if (minVal > cumEnergy[row][minCol + 1]) {
            minIndex = minCol + 1;
        }
        
        return minIndex;
    }
    
    private int minimumVIndex(double[][] cumEnergy, int minRow, int col) {
        int minIndex = minRow -1;
        if (minRow == 0) {
            return cumEnergy[minRow][col] < cumEnergy[minRow + 1][col] ? minRow : minRow + 1;
        }
        
        if (minRow == height - 1) {
            return cumEnergy[minRow - 1][col] < cumEnergy[minRow][col] ? minRow - 1 : minRow;
        }
        
        double minVal = cumEnergy[minRow - 1][col];
        if (minVal > cumEnergy[minRow][col]) {
            minVal = cumEnergy[minRow][col];
            minIndex = minRow;
        }
        
        if (minVal > cumEnergy[minRow + 1][col]) {
            minIndex = minRow + 1;
        }
        
        return minIndex;
    }
    
    private boolean isValidHSeam(int[] seam) {
        for (int j = 0; j < seam.length; ++j) {
            if (seam[j] < 0 || seam[j] > height - 1) {
                return false;
            }
            
            if (j == 0) {
                continue;
            }
            
            if (Math.abs(seam[j] - seam[j - 1]) > 1) {
                return false;
            }
        }
        return true;
    }
    
    private boolean isValidVSeam(int[] seam) {
        for (int j = 0; j < seam.length; ++j) {
            if (seam[j] < 0 || seam[j] > width - 1) {
                return false;
            }
            
            if (j == 0) {
                continue;
            }
            if (Math.abs(seam[j] - seam[j - 1]) > 1) {
                return false;
            }
        }
        return true;
    }
   
   /**
    * current picture
    * @return
    */
   public Picture picture() {
       return picture;
   }
   
   /**
    * width of the current picture
    * @return width
    */
   public int width() {
       return width;
   }
   
   /**
    * height of current picture
    * @return height
    */
   public int height() {
       return height; 
   }
   
   /**
    * energy of pixel at column x and row y
    * @param x
    * @param y
    * @return energy
    */
   public  double energy(int x, int y) {
       double deltaX = 0;
       double deltaY = 0;
       
       if (x < 0 || x > width - 1 || y < 0 || y > height - 1) {
           throw new IllegalArgumentException();
       }
       
       if (0 == x || 0 == y || width - 1 == x || height - 1 == y) {
           return BORDER_ENERGY;
       }
       
       /// Diff for X
       double deltaXR = picture.get(x + 1, y).getRed() 
               - picture.get(x - 1, y).getRed();
       double deltaXG = picture.get(x + 1, y).getGreen() 
               - picture.get(x - 1, y).getGreen();
       double deltaXB = picture.get(x + 1, y).getBlue() 
               - picture.get(x - 1, y).getBlue();
       
       deltaX = Math.pow(deltaXR, 2) 
               + Math.pow(deltaXG, 2) 
               + Math.pow(deltaXB, 2);
       
       /// Diff for Y
       double deltaYR = picture.get(x, y + 1).getRed() 
               - picture.get(x, y - 1).getRed();
       double deltaYG = picture.get(x, y + 1).getGreen() 
               - picture.get(x, y - 1).getGreen();
       double deltaYB = picture.get(x, y + 1).getBlue() 
               - picture.get(x, y - 1).getBlue();
       
       deltaY = Math.pow(deltaYR, 2) 
               + Math.pow(deltaYG, 2) 
               + Math.pow(deltaYB, 2);
       
       return Math.sqrt(deltaX + deltaY);
   }
   
   /**
    * sequence of indices for vertical seam
    * @return
    */
   public   int[] findVerticalSeam() {
       double [][] cumVEnergy = new double[height][width];
       int [] verticalSeam = new int[height];
       
       double minVal = Double.POSITIVE_INFINITY;
       int minIndex = 0;
       
       for (int row = 0; row < height; ++row) {
           for (int col = 0; col < width; ++col) {
               if (0 == row) {
                 cumVEnergy[row][col] = BORDER_ENERGY;
               } else if (0 == col) {
                   cumVEnergy[row][col] = energy(col, row)
                           + Math.min(cumVEnergy[row - 1][col], 
                                   cumVEnergy[row - 1][col + 1]);
               } else if (width - 1 == col) {
                   cumVEnergy[row][col] = energy(col, row)
                           + Math.min(cumVEnergy[row - 1][col - 1], 
                                   cumVEnergy[row - 1][col]);
               } else {
                   cumVEnergy[row][col] = energy(col, row) 
                       + this.min(cumVEnergy[row -1][col - 1], 
                               cumVEnergy[row - 1][col], 
                               cumVEnergy[row - 1][col + 1]);
               }
               if (height - 1 == row) {
                   if (cumVEnergy[row][col] < minVal) {
                       minVal = cumVEnergy[row][col];
                       minIndex = col;
                   }
               }
           }
       }
       
       verticalSeam[height - 1] = minIndex; 
       for (int i = height - 2; i >= 0; --i) {
           verticalSeam[i] = this.minimumHIndex(cumVEnergy, 
                   minIndex,
                   i);
           minIndex = verticalSeam[i];
       }
       
       return verticalSeam;
   }
   
   /**
    * sequence of indices for horizontal seam
    * @return
    */
   public   int[] findHorizontalSeam() {
       double [][] cumHEnergy = new double[height][width];
       int [] horizonSeam = new int[width];
       double minEnergy = Double.POSITIVE_INFINITY;
       int minIndex = 0;
       
       for (int col = 0; col < width; ++col) {
           for (int row = 0; row < height; ++row) {
               if (col == 0) {
                   cumHEnergy[row][col] = BORDER_ENERGY;
               } else if (row == 0) {
                   cumHEnergy[row][col] = energy(col, row) 
                           + Math.min(cumHEnergy[row][col - 1], 
                                   cumHEnergy[row + 1][col - 1]);
               } else if (row == height - 1) {
                   cumHEnergy[row][col] = energy(col, row)
                           + Math.min(cumHEnergy[row - 1][col -1], 
                                   cumHEnergy[row][col -1]);
               } else {
                   cumHEnergy[row][col] = energy(col, row) 
                           + this.min(cumHEnergy[row - 1][col - 1], 
                                   cumHEnergy[row][col - 1], 
                                   cumHEnergy[row + 1][col - 1]);
               }
               
               if (col == width - 1) {
                   if (minEnergy > cumHEnergy[row][col]) {
                       minEnergy = cumHEnergy[row][col];
                       minIndex = row;
                   }
               }
           }
           horizonSeam[width - 1] = minIndex;
           for (int i = width - 2; i >= 0; --i) {
               horizonSeam[i] = this.minimumVIndex(cumHEnergy, minIndex, i);
               minIndex = horizonSeam[i];
           }
       }
       
       return horizonSeam;
   }
   
   /**
    * remove horizontal seam from current picture
    * @param seam
    */
   public void removeHorizontalSeam(int[] seam) {
       if (null == seam) {
           throw new IllegalArgumentException("Horizontal seam is null");
       }
       
       if (seam.length != width()) {
           throw new IllegalArgumentException("Horizontal seam size not matching with picture width");
       }
       
       if (!this.isValidHSeam(seam)) {
           throw new IllegalArgumentException("horizontal Seam Invalid");
       }
       
       Picture horizonCarved = new Picture(width(), height() - 1);
       
       int x = 0;
       for (int col = 0; col < width(); ++col) {
           int y = 0;
           for (int row = 0; row < height(); ++row) {
               if (seam[col] == row) {
                   continue;
               }
               
               horizonCarved.set(x, y, this.picture().get(col, row));
               y++;
           }
           x++;
       }
       
       picture = horizonCarved;
       height = horizonCarved.height();
       width = horizonCarved.width();
   }
   
   /**
    * remove vertical seam from current picture
    * @param seam
    */
   public    void removeVerticalSeam(int[] seam) {
       if (null == seam) {
           throw new IllegalArgumentException("Vertical seam is null");
       }
       
       if (seam.length != height()) {
           throw new IllegalArgumentException("Vertical seam size not matching with picture height");
       }
       
       if (!this.isValidVSeam(seam)) {
           throw new IllegalArgumentException("vertical Seam mismatch");
       }
       
       Picture verticalCarved = new Picture(width() - 1, height());
       int y = 0;
       for (int row = 0; row < height(); ++row) {
           int x = 0;
           for (int col = 0; col < width(); ++col) {
               if (seam[row] == col) {
                   continue;
               }
               
               verticalCarved.set(x, y, this.picture().get(col, row));
               x++;
           }
           y++;
       }
       
       picture = verticalCarved;
       height = picture.height();
       width = picture.width();
   }
   
}
