package sampsonLab;

import java.util.ArrayList;
import java.util.Collections;


/**
 * Created by dfermin on 1/5/17.
 */
public class ESP_Features {

    public double ESP_AA_AC, ESP_EA_AC, ESP_MAF;
    public double ESP_MAX_AA_AE;

    public ESP_Features() {
        ESP_AA_AC = 0;
        ESP_EA_AC = 0;
        ESP_MAF = 0;
        ESP_MAX_AA_AE = 0;
    }


    public static double calcMAF(String AC_str) {
        ArrayList<Double> counts = new ArrayList<Double>();
        ArrayList<Double> AF = new ArrayList<Double>();
        double sum = 0;


        // record all of the allele counts
        for(String s : AC_str.split(",")) {
            double d = Double.valueOf(s);
            counts.add(d);
            sum += d;
        }

        for(double d : counts) AF.add( (d / sum) );

        // Sort AF from high to low and accept the second element as the minor allele frequency (MAF)
        Collections.sort(AF);
        Collections.reverse(AF);

        double MAF_ret = 0;
        if(AF.size() > 1) MAF_ret = AF.get(1);
        else MAF_ret = AF.get(0);

        return MAF_ret;
    }

}
