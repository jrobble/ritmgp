package ritmgp.util;

/**
 * Handles general math functions needed for image processing/matrix math
 *
 * @author Mike Neurohr
 */
public class MathUtil {
    /**
     * standard deviation of the normal distribution used by the MathCAD
     * ksmooth function
     * http://www.imakenews.com/ptcexpress/e_article001119634.cfm?x=bcR12Vy,b3jsqcsB,w
     */
    private static final double STD_DEV = 0.3708158988058244;

    /**
     * Calulates the plain average of an array
     *
     * @param vec array to average
     * @return the average of vec
     */
    public static double avg(double[] vec) {
        if (vec.length == 0) { //if there are no elements, the average is 0
            return 0;
        }
        double sum = 0;
        for (double elem : vec) { //for each element in the vector
            sum += elem; //accumulate them
        }
        return sum / vec.length;
    }

    /**
     * Calulates a linear interpolation of a third point between two others.
     *
     * @param x2 x value of the first point
     * @param y2 y value of the first point
     * @param x1 x value of the second point
     * @param y1 y value of the second point
     * @param x x value to interpolate a y value for
     * @return the y value at x based on (x1, y1) and (x2, y2)
     */
    private static double linterp(double x2, double y2, double x1, double y1, double x) {
        double m = (y2 - y1) / (x2 - x1);
        double b = y1 - m * x1;
        return m * x + b;
    }

    public static double linterp(double[] xs, double[] ys, double x){
        int index = 0;
        for(;index < xs.length - 1; index++){
            if(xs[index] <= x && xs[index + 1] > x) break;
        }
        if(xs[index] == x) return ys[index]; //exact match
        if(index == xs.length - 1) index--;
        return linterp(xs[index], ys[index], xs[index + 1], ys[index + 1], x);
    }

    /**
     * Does a kernal smoothing on a given array with a Gaussian kernel
     * @param array array to smooth
     * @param bw bandwidth to use for Gaussian smoothing
     * @return an array of equal size which is smoothed via Gaussian kernel
     */
    public static double[] ksmooth(double[] array, double bw) {
        double[] ans = new double[array.length]; //placeholder for the answer
        for (int i = 0; i < array.length; i++) { //we need to fill in each element in the array
            double num = 0; //numerator
            double denom = 0; //denominator

            //do both summations at the same time
            for (int j = 0; j < array.length; j++) {
                num += getGaussWeight((i - j) / bw) * array[j];
                denom += getGaussWeight((i - j) / bw);
            }

            //final division after summations
            ans[i] = num / denom;
        }
        return ans;
    }

    /**
     * Gives the value on a bell curve for a given x
     * @param x value to plug into the bell curve formula
     * @return the y value on the bell curve at x
     */
    private static double getGaussWeight(double x) {
        return 1 / ((2 * Math.PI) * STD_DEV) * Math.exp(-x * x / (2 * STD_DEV * STD_DEV));
    }

    public static double max(double[] array){
        double max = array[0];
        for(double elem : array){
            if(elem > max) max = elem;
        }
        
        return max;
    }

    public static double min(double[] array){
        double min = array[0];
        for(double elem : array){
            if(elem < min) min = elem;
        }

        return min;
    }

    public static int min(int[] array){
        int min = array[0];
        for(int elem : array){
            if(elem < min) min = elem;
        }
        return min;
    }

    public static int max(int[] array){
        int max = array[0];
        for(int elem : array){
            if(elem > max) max = elem;
        }

        return max;
    }
}
