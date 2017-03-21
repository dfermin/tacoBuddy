package sampsonLab;

/**
 * Created by dfermin on 3/15/17.
 */

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;

public class VCFDB {
    public static Connection conn = null;
    public static String dbName;

    public VCFDB(String s) throws SQLException {

        String URL = "jdbc:hsqldb:mem:vcfdb";
        if((null != s) && (s.length() > 1) ) {
            dbName = s;
            URL = "jdbc:hsqldb:hsql://localhost/" + s;
        }
        else dbName = s;

        try {
            // Registering the HSQLDB JDBC driver
            Class.forName("org.hsqldb.jdbc.JDBCDriver");

            // Create connection with HSQLDB
            conn = DriverManager.getConnection(URL, "SA", "");
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
    }


    /**********************************************************************************/
    // Function creates a table just for the features we are interested in
    public void createVCFtable() {

        String query = "CREATE TABLE vcf (" +
                "id VARCHAR(100) NOT NULL";  // this is the variant coordinate: 'chr':'pos'

        // Add the features the user asked for that are stored in 'globals.featureMap'
//        for(String k : globalFunctions.featureSet) {
//
//        }
    }

}
