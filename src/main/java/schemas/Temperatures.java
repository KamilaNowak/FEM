package schemas;

import lombok.*;

@Getter
@Setter
@AllArgsConstructor
@NoArgsConstructor
@ToString
public class Temperatures {

    double time;
    double minTemperature;
    double maxTemperature;
}
