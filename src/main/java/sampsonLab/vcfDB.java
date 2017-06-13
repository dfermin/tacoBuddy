package sampsonLab;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.Statement;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Created by dfermin on 6/7/17.
 */
public class vcfDB {

    Connection conn = null;
    SortedSet<String> columns = null;

    public vcfDB(boolean inMemory) {

        try{
            Class.forName("com.mysql.jdbc.Driver");
            conn = DriverManager.getConnection("jdbc:mysql://localhost/dfermin?user=dfermin&password=hsp101");
        } catch (SQLException e) {
            e.printStackTrace();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }

        columns = new TreeSet<String>();
    }


    // Function to create the initial tables
    public void createVCF(SortedSet<String> featureSet) throws SQLException {

        /*##############################################################################################################
        ** We create 1 table per chromosome. The table may be empty if no variants are found for it in the VCF file
        **############################################################################################################*/
        Statement stmt = conn.createStatement();

        String q = "CREATE TABLE variants ("
                 + "  chr VARCHAR(2) NOT NULL, "
                 + "  pos INT NOT NULL, "
                 + "  allowedSite TINYINT, "
                 + "  geneId VARCHAR(250), "
                 + "  dbSNPID TEXT, "
                 + "  ref TEXT, "
                 + "  alt TEXT, "
                 + "  sampleMAF DOUBLE, ";


        // Append the user's requested features to the primary table string
        for(String s : featureSet) {

            if(s.equalsIgnoreCase("sample_maf")) continue;

            columns.add(s.toUpperCase());

            // This case is special and must be manually defined
            if(s.equalsIgnoreCase("EFF")) {
                q += "  EFF TEXT,  ";
            }

            if(
                s.equalsIgnoreCase("MUTATIONTASTER") ||
                s.equalsIgnoreCase("POLYPHEN2_HVAR") ||
                s.equalsIgnoreCase("SIFT") ||
                s.equalsIgnoreCase("GERP")) {

                dbNSFP_Features tmp = new dbNSFP_Features();
                q += s + " " + tmp.getDataType(s) + ",  ";
            }

            if(s.toUpperCase().startsWith("ESP_")) {
                ESP_Features tmp = new ESP_Features();
                q += s + " " + tmp.getDataType(s) + ",  ";
            }

            if(s.equalsIgnoreCase("IS_LOF")) {
                q += s + " TINYINT,  ";
            }
        }
        String tmp = trimEnd(q, 3);
        q = tmp + " ); ";
        tmp = null;


        // Iterate over the chromosomes and create a table foreach
        //String[] genome = { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "Y" };

        stmt.executeUpdate("DROP TABLE IF EXISTS variants; ");
        stmt.executeUpdate(q);
        stmt.executeUpdate("ALTER TABLE variants ADD INDEX (chr, pos)");
        stmt.executeUpdate("ALTER TABLE variants ADD INDEX (geneID)");
    }


    // Function takes a passed VariantInfo object and loads it's featureSet data into the database
    public void loadVariant(VariantInfo vi, String geneId) throws SQLException {


        String load_stmt = "INSERT INTO variants (chr, pos, allowedSite, geneId, dbSNPID, ref, alt, sampleMAF, ";

        for(String s : columns) { load_stmt += s + ", "; }
        load_stmt = trimEnd(load_stmt, 2); // removes the trailing ", ";

        load_stmt += ") VALUES (";
        load_stmt += " '" + vi.chr + "', " +
                String.valueOf(vi.pos) + ", " +
                vi.allowedSite + ", " +
                "'" + geneId + "', " +
                "'" + vi.dbsnp_id + "', " +
                "'" + vi.REF + "', " +
                "'" + vi.ALT + "', " +
                String.valueOf(vi.sample_MAF) + ", ";

        for(String s : columns) {
            load_stmt += vi.getDBformattedFeature(s) + ", ";
        }
        load_stmt = trimEnd(load_stmt, 2);
        load_stmt += "); ";

        conn.createStatement().executeUpdate(load_stmt);


    }

    /*************************************************************************************************/
    public static String trimEnd(String s, int N) {
        // s = string to trim
        // N = number of characters to remove from the end

        int len = s.length() - N;
        return s.substring(0, len);
    }
}
