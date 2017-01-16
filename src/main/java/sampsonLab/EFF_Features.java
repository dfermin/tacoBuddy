package sampsonLab;

import com.google.common.base.Joiner;

import java.util.*;

/**
 * Created by dfermin on 1/6/17.
 */
public class EFF_Features extends FeatureClass {

    public HashMap<String, String> eff; // k = transcript ID, v = status_for_this_transcript

    public EFF_Features(String ss) {

        eff = new HashMap<String, String>();

        String[] inputAry = ss.split(",");

        for(String s : inputAry) {

            String[] tmpAry = s.split("\\|");

            String label = tmpAry[0].substring(0,tmpAry[0].indexOf('(')).trim(); // gets the "assigned effect" of this variant on the given gene.
            String transId = "";
            for(int j = 1; j < tmpAry.length; j++) {
                if( tmpAry[j].startsWith("ENST00") ) {
                    transId = tmpAry[j];
                    break;
                }
            }
            if(transId.startsWith("ENST00")) eff.put(transId, label);
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


    // Function iterates over the values in 'eff' and concatenates the unique ones into an array
    // that jexl can interpret. This function is for gene-level output.
    public String returnJexlArray() {
        String ret = "";
        SortedSet<String> SS = new TreeSet<String>();

        for(String v : eff.values()) {

            for(String i : v.split(",")) {
                i = "'" + i.trim() + "'";
                SS.add(i);
            }
        }

        ret = "[" + Joiner.on(",").join(SS) + "]";
        return ret;
    }
}
