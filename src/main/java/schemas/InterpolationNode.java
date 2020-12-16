package schemas;

import lombok.*;

@Getter
@Setter
@AllArgsConstructor
@NoArgsConstructor
@ToString

public class InterpolationNode {
    private double ksi;
    private double eta;
    private double weight;

    public InterpolationNode(double ksi, double eta) {
        this.ksi = ksi;
        this.eta = eta;
    }
}
