package schemas;

public class InterpolationNode {
    private double ksi;
    private double eta;
    private double weight;

    public InterpolationNode(double ksi, double eta) {
        this.ksi = ksi;
        this.eta = eta;
    }

    public InterpolationNode() {
    }

    public double getKsi() {
        return ksi;
    }

    public void setKsi(double ksi) {
        this.ksi = ksi;
    }

    public double getEta() {
        return eta;
    }

    public void setEta(double eta) {
        this.eta = eta;
    }

    public double getWeight() {
        return weight;
    }

    public void setWeight(double weight) {
        this.weight = weight;
    }

    @Override
    public String toString() {
        return "InterpolationNode{" +
                "ksi=" + ksi +
                ", eta=" + eta +
                ", weight=" + weight +
                '}';
    }
}
