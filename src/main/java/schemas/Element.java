package schemas;

import lombok.*;
import java.util.List;

@Getter
@Setter
@AllArgsConstructor
@NoArgsConstructor
@ToString
public class Element {

    private List<Integer> IDs;
    private List<Node> nodes;
    private double[][] H = new double[4][4];
    private double[][] C = new double[4][4];
    private double[][] Hbc = new double[4][4];
    private double[] P = new double[4];
}
