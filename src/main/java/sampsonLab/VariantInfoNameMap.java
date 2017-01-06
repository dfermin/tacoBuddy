package sampsonLab;

import java.util.HashMap;

/**
 * Created by dfermin on 1/6/17.
 *
 * This file maps between the "observed" or "reported" names for the variant features and their
 * "common" names. For instance GERP++NR = "dbNSFP_GERP__NR"
 *
 * Syntax: key = "common", value = "reported"
 */

public class VariantInfoNameMap {
    public static HashMap<String, String> fieldMap = null;

    public VariantInfoNameMap() {
        fieldMap = new HashMap<String, String>();

        fieldMap.put("GERP++".toUpperCase(),       "dbNSFP_GERP___RS");
        fieldMap.put("POLYPHEN2_HVAR".toUpperCase(), "dbNSFP_Polyphen2_HVAR_pred");
        fieldMap.put("SIFT".toUpperCase(),           "dbNSFP_SIFT_pred");
        fieldMap.put("MUTATIONTASTER".toUpperCase(), "dbNSFP_MutationTaster_pred");
        fieldMap.put("is_lof".toUpperCase(),         "EFF");
    }


    public static boolean inMap(String k) {
        return fieldMap.containsKey(k);
    }

    public static String getValue(String k) {

        if(fieldMap.containsKey(k)) {
            return fieldMap.get(k);
        }
        return null;
    }
}
