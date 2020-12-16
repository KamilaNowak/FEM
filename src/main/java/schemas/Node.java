package schemas;

import lombok.*;

@Getter
@Setter
@AllArgsConstructor
@NoArgsConstructor
@ToString
public class Node {
    private int id;
    private double x;
    private double y;
    private double t;
    private boolean bc = false;

    public Node(double x, double y) {
        this.x = x;
        this.y = y;
    }
}
