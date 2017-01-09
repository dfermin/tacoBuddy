package sampsonLab;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

/**
 * Created by dfermin on 1/6/17.
 */
public class EFF_Features extends FeatureClass {

    public HashMap<String, String> eff; // k = transcript ID, v = status_for_this_transcript

    public EFF_Features(Object o) {

        eff = new HashMap<String, String>();

        ArrayList<String> inputAry = (ArrayList<String>) o;

        for(String s : inputAry) {

            String[] tmpAry = s.split("\\|");

            String label = tmpAry[0].substring(0,tmpAry[0].indexOf('(')); // gets the "assigned effect" of this variant on the given gene.
            String transId = "";
            for(int j = 1; j < tmpAry.length; j++) {
                if( tmpAry[j].startsWith("ENST00") ) {
                    transId = tmpAry[j];
                    break;
                }
            }
            eff.put(transId, label);
        }
    }


    // Function returns the value of 'eff' map for the given key 'k'.
    // If no match is found, null is returned
    public String findTS(String search_str) {
        String ret = "no_match";

        if(eff.containsKey(search_str)) {
            ret = eff.get(search_str);
        }
        return ret;
    }

}
