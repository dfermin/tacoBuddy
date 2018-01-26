package sampsonLab;

import com.sampsonlab.filter.JavascriptFilter;

import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;
import java.util.HashMap;


/**
 * Created by dfermin on 1/25/18.
 *
 * This class uses the Nashorn Javascript engine built into JAVA 8 to evaluate the filtering expressions
 * given by the user
 */
public class VariantFiltering {
    JavascriptFilter filter; // Actual object that performs the filtering
    String filterStr = null; // The filter string you want evaluated

    public VariantFiltering(String filterPattern) {
        ScriptEngine engine = new ScriptEngineManager().getEngineByName("nashorn");
        filterStr = filterPattern.replaceAll(" and ", " && ");
        filterStr = filterStr.replace(" or ", " || ");

        filter = JavascriptFilter.createWithEngine(engine);
        filter.setDefaultFilterStr(filterStr);
    }


    public boolean evalVariant(HashMap<String, Object> params) throws ScriptException {
        boolean ret = filter.apply(params);
        return ret;
    }

}
