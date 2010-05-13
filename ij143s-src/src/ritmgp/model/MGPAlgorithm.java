package ritmgp.model;

import ij.ImagePlus;
import ritmgp.util.MathUtil;
import ritmgp.util.Matrix;

/**
 *
 * @author Mike Neurohr
 */
public class MGPAlgorithm {
    private ImagePlus bright;
    private ImagePlus dark;
    private double refArea;
    private double refReflect;
    private double fov;     //fov and exposure will need to be worked out
                            //on the device, they may need software conversions
    private double fStop;
    private boolean autoBaseline;
    private double leftAngle;
    private double rightAngle;

    private static double workingDist = 183.27; //183.27 for our device
    private static double instAngle = 18.25; //19 degrees for our device?
    private static double cylDiam = 30; //cylinder diameter

    public MGPAlgorithm(ImagePlus bright, ImagePlus dark, double refArea,
            double refReflect, double fov, double fStop, boolean autoBaseline,
            double leftAngle, double rightAngle) {
        this.bright = bright;
        this.dark = dark;
        this.fStop = fStop;
        this.refArea = refArea;
        this.refReflect = refReflect;
        this.autoBaseline = autoBaseline;
        this.leftAngle = leftAngle;
        this.rightAngle = rightAngle;
        this.fov = fov;
    }
    
    /**
     * The algorithm coded here is copied directly from the MathCAD code
     * in the docs folder of this project. Changes are noted inline with code.
     * 
     * @return A result object holding the result vector, the baseline angles
     * used, the final BRDF, and the angles corresponding to the BRDF
     */
    public MGPResult measureGloss(){
        Matrix brightMat = new Matrix(getPixels(bright));
        Matrix darkMat = new Matrix(getPixels(dark));

        if(brightMat.getCols() != darkMat.getCols() ||
                brightMat.getRows() != darkMat.getRows())
            return null;

        int numCols = brightMat.getCols();
        int numRows = brightMat.getRows();
        //irradiance difference
        Matrix diff0 = Matrix.subtract(brightMat, darkMat); //deltaI0 in MathCAD
        diff0 = Matrix.divideScalar(diff0, fStop * 1); //gain is 1

        Matrix diff = Matrix.zeroNegatives(diff0);
        //irradiance difference
        
        //correct for illumination variations
        double[] sq = new double[numRows];

        for(int i = 0; i < sq.length; i++){
            sq[i] = MathUtil.avg(diff.selectRow(i));
        }
        
        double[] ff = MathUtil.ksmooth(sq, diff.getRows() / 5);
        
        double meanFF = MathUtil.avg(ff);

        //corresponds to deltaIF in MathCAD code
        double[][] diffFArr = new double[numRows][numCols];

        for(int i = 0; i < diffFArr.length; i++){
            for(int j = 0; j < diffFArr[0].length;j++){
                diffFArr[i][j] = diff.get(i, j) / ff[i] * meanFF;
            }
        }

        Matrix diffF = new Matrix(diffFArr);
        //SqF was never used after its creation, so it is not coded here

        //correct for illumination variations
        
        //raw BRDF and noise functions
        double[] BRDF0 = new double[numCols];

        for(int j = 0; j < BRDF0.length; j++){
            BRDF0[j] = MathUtil.avg(diffF.selectCol(j));
        }
        
        //meanBRDF0 and maxBRDF0 were not used so they are not calculated here
        double minBRDF0 = MathUtil.min(BRDF0);

        double[] BRDF1 = new double[numCols];

        for(int j = 0; j < BRDF0.length; j++){
            BRDF1[j] = BRDF0[j] - minBRDF0;
        }

        double maxBRDF1 = MathUtil.max(BRDF1);

        double[] noise1 = new double[BRDF1.length];

        for(int j = 0; j < noise1.length; j++){
            double sum = 0;
            for(int i = 0; i < diffF.getRows();i++){
                sum += Math.pow(diffF.get(i, j) - BRDF0[j], 2);
            }
            noise1[j] = Math.sqrt(sum / diffF.getRows()) / maxBRDF1;
        }
        int jp = numCols / 2;
        //raw BRDF and noise functions

        //calculate the angles
        double errorAngle  = Math.atan(fov/(workingDist * 2)) * 180 / Math.PI;
        double k = (numCols - 1) / fov;
        final double ap = 1000 / k;//part of final result
        double fresnells = instAngle / 2;//corresponds to theta in MathCAD
        double jc = (k * fov) / 2 * Math.sin(fresnells * Math.PI / 180) + jp;
        
        double[] beta1 = new double[numCols];
        double[] beta = new double[numCols];
        double[] alpha1 = new double[numCols];
        for(int j = 0; j < numCols; j++){
            beta1[j] = Math.asin(2 * (jc - j) / (k * cylDiam)) * 180 / Math.PI;
            beta[j] = beta1[j] * 1 / ((errorAngle / 90) + 1);
            alpha1[j] = fresnells - beta[j];
        }
        //calculate the angles

        //eliminate imaginary angles
        int[] IM = new int[numCols];

        for(int j = 0; j < IM.length; j++){
            IM[j] =  Double.isNaN(alpha1[j]) ? numCols - 1 : j;
        }
        
        int JIM = MathUtil.min(IM);

        double[] alpha2 = new double[alpha1.length];

        for(int j = 0; j < alpha2.length; j++){
            alpha2[j] = (j < JIM) ? alpha1[JIM] : alpha1[j];
        }
        
        double[] BRDF2 = new double[numCols];

        for(int j = 0; j < BRDF2.length; j++){
            BRDF2[j] = (j < JIM) ? BRDF1[JIM] : BRDF1[j];
        }
        
        double[] noise2 = new double[numCols];

        for(int j = 0; j < noise2.length; j++){
            noise2[j] = (j < JIM) ? noise1[JIM] : noise1[j];
        }
        
        //MaxBRDF2 was not used in MathCAD so it is not coded here
        //eliminate imaginary angles

        //Generate functions at equally spaced angles
        //alpha3 was created for MathCAD trickery
        //it is not used here
        double alpha2Min = MathUtil.min(alpha2);
        double alpha2Max = MathUtil.max(alpha2);

        double alphaInt = (alpha2Max - alpha2Min) / numCols;

        double[] alpha = new double[numCols];
        double[] BRDF4 = new double[numCols]; //BRDF3 skipped because it is
                                              //unnecessary in this program but
                                              //names should be consistent
         double[] noise = new double[numCols];//noise3 skipped because it is
                                              //unnecessary in this program but
                                              //names should be consistent
        for(int j = 0; j < numCols; j++){
            alpha[j] = j * alphaInt + alpha2Min;
            BRDF4[j] = MathUtil.linterp(alpha2, BRDF2, alpha[j]);
            noise[j] = MathUtil.linterp(alpha2, noise2, alpha[j]);
        }

        double alphaMin = MathUtil.min(alpha);
        double alphaMax = MathUtil.max(alpha);
        double maxBRDF4 = MathUtil.max(BRDF4);
        double maxBRDF420 = maxBRDF4 / 20;
//        int maxBRDF4j = 0;
//
//        for(int j = 0; j < numCols; j++){
//            if(BRDF4[j] == maxBRDF4){
//                maxBRDF4j = j;
//            }
//        }
//
//        double maxBRDF4alpha = alpha[maxBRDF4j];
        final double granularity = MathUtil.max(noise);
        //Generate functions at equally spaced angles

        //auto baseline
        double alphaLeft, alphaRight;

        int[] JLV = new int[numCols];
        for (int j = 0; j < numCols; j++) {
            if (alpha[j] > 0) {
                JLV[j] = 0;
            } else if (BRDF4[j] > maxBRDF420) {
                JLV[j] = 0;
            } else {
                JLV[j] = j;
            }
        }
        int JLv = MathUtil.max(JLV);

        double alphaX = 2 * (1 + Math.exp(-0.5 * Math.abs(alpha[JLv])))
                * Math.abs(alpha[JLv]);

        double alphaLimit = Math.abs(alphaMin) < alphaMax ? Math.abs(alphaMin) : alphaMax;

        double alphaBase = alphaX > alphaLimit ? alphaLimit : alphaX;

        int[] JWL = new int[numCols];
        int[] JWR = new int[numCols];
        for (int j = 0; j < numCols; j++) {
            if (alpha[j] > 0) {
                JWL[j] = 0;
            } else if (alpha[j] > -alphaBase) {
                JWL[j] = 0;
            } else {
                JWL[j] = j;
            }
        }
        
        int JL = Math.min(MathUtil.max(JWL) + 2, numCols);
        for (int j = 0; j < numCols; j++) {
            if (alpha[j] < 0) {
                JWR[j] = 0;
            } else if (alpha[j] > -alpha[JL]) {
                JWR[j] = 0;
            } else {
                JWR[j] = j;
            }
        }
        int JR = Math.max(MathUtil.max(JWR) - 2, 0);
        //auto baseline

        //baseline correction
        if (autoBaseline) {
            alphaLeft = alpha[JL];
            alphaRight = alpha[JR];
        } else {
            alphaLeft = leftAngle;
            alphaRight = rightAngle;
        }

        double vLeft = BRDF4[JL];
        double vRight = BRDF4[JR];

        double[] X = {alphaLeft, alphaRight};
        double[] Y = {vLeft, vRight};

        double[] VB = new double[numCols];
        double[] BRDF5 = new double[numCols];
        double[] BRDF6 = new double[numCols];
        double[] BRDF77 = new double[numCols];
        double[] BRDF7 = new double[numCols];
        for(int j = 0; j < numCols; j++){
            VB[j] = MathUtil.linterp(X, Y, alpha[j]);
            BRDF5[j] = BRDF4[j] - VB[j];
            if (alpha[j] < alphaLeft) {
                BRDF6[j] = 0;
            } else if (alpha[j] > alphaRight) {
                BRDF6[j] = 0;
            } else {
                BRDF6[j] = BRDF5[j];
            }
            BRDF77[j] = BRDF6[j] - VB[j];
            BRDF7[j] = BRDF77[j] < 0 ? 0 : BRDF77[j];
        }
        //baseline correction
        
        //calibrating the BRDF
        double A0 = 0;
        double[] BRDF8 = new double[numCols];
        for(int j = 0;j < BRDF7.length; j++){
            A0 += BRDF7[j];
            BRDF8[j] = BRDF7[j] / refArea;
        }

        final double A = A0 * alphaInt / refArea;
        //calibrating the BRDF

        //stats
        final double h = MathUtil.max(BRDF8);
        final double rho = A * refReflect;

        double hHalf = h/2;

        int[] halfLV = new int[numCols];
        int[] halfRV = new int[numCols];

        double h10 = h / 10;

        int[] tenthLV = new int[numCols];
        int[] tenthRV = new int[numCols];
        for(int j = 0;j < halfLV.length; j++){
            if(alpha[j] > 0){
                halfLV[j] = 0;
            }else if(BRDF8[j] > hHalf){
                halfLV[j] = 0;
            }else{
                halfLV[j] = j;
            }

            if(alpha[j] < 0){
                halfRV[j] = 0;
            }else if(BRDF8[j] > hHalf){
                halfRV[j] = j;
            }else{
                halfRV[j] = 0;
            }

             if(alpha[j] > 0){
                tenthLV[j] = 0;
            }else if(BRDF8[j] > h10){
                tenthLV[j] = 0;
            }else{
                tenthLV[j] = j;
            }

            if(alpha[j] < 0){
                tenthRV[j] = 0;
            }else if(BRDF8[j] > h10){
                tenthRV[j] = j;
            }else{
                tenthRV[j] = 0;
            }
        }

        int halfL = MathUtil.max(halfLV);
        int halfR = MathUtil.max(halfRV);

        //VwLeft and VwRight were only used to create Yw which was not used
        double alphaWL = alpha[halfL];
        double alphaWR = alpha[halfR];

        //Xw and Yw were not used
        final double wHalf = alphaWR - alphaWL;

        int tenthL = MathUtil.max(tenthLV);
        int tenthR = MathUtil.max(tenthRV);

        //v10Left and v10Right were only used to create Y10 which was not used
        double alpha10L = alpha[tenthL];
        double alpha10R = alpha[tenthR];

        //X10 and Y10 were not used
        final double w10 = alpha10R - alpha10L;

        double M2 = 0, M3 = 0, M4 = 0;
        double[] PDF = new double[numCols];
        for(int j = 0; j < PDF.length; j++){
            PDF[j] = BRDF8[j] / A;
            M2 += PDF[j] * alpha[j] * alpha[j];
            M3 += PDF[j] * Math.pow(alpha[j], 3);
            M4 += PDF[j] * Math.pow(alpha[j], 4);
        }
        final double sigma = Math.sqrt(M2);
        final double skewness = M3 / Math.pow(M2, 3.0 / 2);
        final double kurtosis = M4 / (M2 * M2) - 3;
        //stats

        Double[] results = {A, wHalf, w10, h, rho, granularity, ap, sigma,
                            skewness, kurtosis};
        return new MGPResult(results, BRDF4, alphaLeft, alphaRight, alpha, BRDF8);
    }
    /**
     * Converts an image from an ImagePlus to a double[][] of
     * grayscale pixel values.
     *
     * @param image the image to convert
     * @return a 2D array of grayscale pixel values
     */
    private double[][] getPixels(ImagePlus image) {
        double[][] ans = new double[image.getHeight()][image.getWidth()];
        for(int i = 0; i < image.getHeight(); i++){
            for(int j = 0; j < image.getWidth(); j++){
                ans[i][j] = image.getPixel(j, i)[0]; //getPixel takes x, y
                                                    //which is backwards from
                                                    //row, col
            }
        }
        return ans;
    }
}
