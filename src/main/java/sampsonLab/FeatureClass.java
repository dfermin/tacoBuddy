package sampsonLab;

import java.util.HashMap;

/**
 * Created by dfermin on 1/9/17.
 */

// This is the generic class. The specific VCF Info features are inherited from this class
public class FeatureClass {

    // Theoretically, variables common to all the the VCF Info features would go here.
    // But as fo 2016.Jan.09 all of our Features are so different that we have nothing in common to put here.
    // Honestly this class was created so I could dynamically store the various feature objects.

    public HashMap<String, String> dataTypeMap = null; // holds SQL mapping for data types.

    public void setDataType(String k, String v) {
        if(null == dataTypeMap) dataTypeMap = new HashMap<String, String>();
        dataTypeMap.put(k, v);
    }

    public String getDataType(String k) {
        String ret = "#NULL";

        if(dataTypeMap.containsKey(k)) {
            ret = dataTypeMap.get(k);
        }
        else {
            System.err.println("\nERROR! Unable to find '" + k + "' in the FeatureClass::dataTypeMap object for " + this.getClass().getSimpleName());
            System.exit(-1);
        }
        return(ret);
    }
}
