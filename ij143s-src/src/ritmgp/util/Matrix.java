package ritmgp.util;

/**
 *A wrapper class for a 2-D array to handle simple matrix math operations
 *
 * @author Mike Neurohr
 */
public class Matrix {
    private double[][] mat; //the actual array of numbers

    public Matrix(double[][] mat){
        this.mat = mat;
    }

    /**
     * @return the number of rows in this matrix
     */
    public int getRows(){
        return mat.length;
    }

    /**
     * @return the number of columns in this matrix
     */
    public int getCols(){
        return mat[0].length;
    }

    public double get(int i, int j){
        return mat[i][j];
    }

    /**
     * Select a row by index from this matrix
     *
     * @param index the index of the row to select
     * @return the row of the array
     */
    public double[] selectRow(int index){
        return mat[index];
    }

    /**
     * Select a column by index from this matrix
     *
     * @param index the index of the row to select
     * @return the column of the array
     */
    public double[] selectCol(int index){
        double[] ans = new double[getRows()];
        for(int i = 0; i < ans.length; i++){
            ans[i] = mat[i][index];
        }
        return ans;
    }

    /**
     * Element by element subtraction of two matrices
     *
     * @param m1 the first matrix
     * @param m2 the second matrix
     * @return a new matrix which is the result of an element-wise subraction
     * of m1 and m2
     */
    public static Matrix subtract(Matrix m1, Matrix m2){
        if(m1.getRows() != m2.getRows() || m1.getCols() != m2.getCols()) return null;
        double[][] ans = new double[m1.getRows()][m1.getCols()];
        for(int row = 0; row < ans.length;row++){
            for(int col = 0; col < ans[0].length;col++){
                ans[row][col] = m1.mat[row][col] - m2.mat[row][col];
            }
        }
        return new Matrix(ans);
    }

    /**
     * Scalar division of a matrix
     *
     * @param m the matrix to divide
     * @param scalar the scalar to divide by
     * @return a new matrix which is the result of scalar division of
     * m by scalar
     */
    public static Matrix divideScalar(Matrix m, double scalar){
        double[][] ans = new double[m.getRows()][m.getCols()];
        for(int row = 0; row < ans.length;row++){
            for(int col = 0; col < ans[0].length;col++){
                ans[row][col] = m.mat[row][col] / scalar;
            }
        }
        return new Matrix(ans);
    }

    /**
     * Sets negative values in a Matrix to 0
     *
     * @param m the matrix to change
     * @return a new matrix which is a copy of m, but with negative
     * values set to 0.
     */
    public static Matrix zeroNegatives(Matrix m){
     double[][] ans = new double[m.getRows()][m.getCols()];
        for(int row = 0; row < ans.length;row++){
            for(int col = 0; col < ans[0].length;col++){
                ans[row][col] = (m.mat[row][col] < 0) ? 0 : m.mat[row][col];
            }
        }
        return new Matrix(ans);
    }
}
