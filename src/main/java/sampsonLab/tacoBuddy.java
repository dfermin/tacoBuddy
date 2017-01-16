/**
 * Created by dfermin on 11/29/16.
 */

package sampsonLab;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Map;

import apple.laf.JRSUIConstants;
import com.google.common.base.Joiner;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.apache.commons.cli.ParseException;


public class tacoBuddy {


    static public globalFunctions globals;
    static public VariantInfoNameMap VCF_Info_Name_Map;
    static VCFParser the_vcf_parser;
    static ArrayList<JRSUIConstants.Variant> DOM_var = null, REC_var = null;


    public static void main(String[] args) throws IOException, ParseException, SQLException {

        // Prepare the globals object to handle user input
        globals = new globalFunctions();
        VCF_Info_Name_Map = new VariantInfoNameMap();


        if( args.length < 1 ) {
            System.err.print("\nUSAGE: java -jar tacoBuddy.jar -i <input_file> or -t (to generate template input file)\n\n");
            return;
        }

        globals.parseCommandLineArgs(args);
        globals.parseGFF3();

        the_vcf_parser = new VCFParser(globals.inputVCF, globals.inputVCF_tabix);


        // Print header line
        String hdr = "#SAMPLE\tGenotype\tFilter\tGene_Name\tCoord (hg19)\t";
        hdr += Joiner.on("\t").join(globalFunctions.featureSet);
        hdr = hdr.replaceAll("EFF\t", "");
        System.out.println(hdr);

        // Record any variants that land within the exons of the genes stored in DOM and REC maps
        if(globals.DOM_geneMap.size() > 0) the_vcf_parser.parse(globals.DOM_geneMap, globals.filterDOM, "DOM");
        if(globals.REC_geneMap.size() > 0) the_vcf_parser.parse(globals.REC_geneMap, globals.filterREC, "REC");

    }

}
