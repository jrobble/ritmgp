package ritmgp.model;

/**
 * A neat little package for the measurment result vector and the final BRDF
 *
 * @author Mike Neurohr
 */
public class MGPResult {
    private Double[] resVector;
    private double[] BRDF;
    private double leftAngle;
    private double rightAngle;
    private double[] alpha;

    MGPResult(Double[] results, double[] BRDF, double leftAngle,
            double rightAngle, double[] alpha) {
        resVector = results;
        this.BRDF = BRDF;
        this.leftAngle = leftAngle;
        this.rightAngle = rightAngle;
        this.alpha = alpha;
    }

    public Double[] getVector() {
        return resVector;
    }

    public double[] getBRDF() {
        return BRDF;
    }

    public double getLeftAngle() {
        return leftAngle;
    }

    public double getRightAngle() {
        return rightAngle;
    }

    public double[] getAlpha(){
        return alpha;
    }
}
