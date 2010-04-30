package ritmgp.model;

/**
 * A neat little package for the measurment result vector and the final BRDF
 *
 * @author Mike Neurohr
 */
public class MGPResult {
    private Double[] resVector;
    private double[] graphBRDF;
    private double[] finalBRDF;
    private double leftAngle;
    private double rightAngle;
    private double[] alpha;

    MGPResult(Double[] results, double[] graphBRDF, double leftAngle,
            double rightAngle, double[] alpha, double[] finalBRDF) {
        resVector = results;
        this.graphBRDF = graphBRDF;
        this.finalBRDF = finalBRDF;
        this.leftAngle = leftAngle;
        this.rightAngle = rightAngle;
        this.alpha = alpha;
    }

    public Double[] getVector() {
        return resVector;
    }

    public double[] getGraphBRDF() {
        return graphBRDF;
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

    public double[] getFinalBRDF() {
        return finalBRDF;
    }
}
