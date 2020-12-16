package schemas;

import lombok.*;
import java.util.ArrayList;
import java.util.List;

@Getter
@Setter
@AllArgsConstructor
@NoArgsConstructor
@ToString
public class Network {

    private List<Node> nodes = new ArrayList<>();
    private List<Element> elements = new ArrayList<>();

    public void printNodes() {
        for ( int i =0;i< nodes.size(); i++)
            System.out.println(i+" : (" + nodes.get(i).getX() + ", " + nodes.get(i).getY() + ")");
    }

    public void printElements() {
        for (Element element : elements) System.out.print("\n[" + element + "]");
    }
}
