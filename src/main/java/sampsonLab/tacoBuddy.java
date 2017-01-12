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
        String hdr = "passedFilter\tCoord (hg19)\tgeneid\t";
        hdr += Joiner.on("\t").join(globalFunctions.featureSet);
        System.out.println(hdr);

        // Record any variants that land within the exons of the genes stored in DOM and REC maps
        if(globals.DOM_geneMap.size() > 0) the_vcf_parser.parse(globals.DOM_geneMap, globals.filterDOM);
        if(globals.REC_geneMap.size() > 0) the_vcf_parser.parse(globals.REC_geneMap, globals.filterREC);

        //the_vcf_parser.db_variants();


        System.exit(0);


        // This command only works if you have the tabix tbi file for the input VCF
        VCFFileReader vcfr = new VCFFileReader(globals.inputVCF, globals.inputVCF_tabix);

        CloseableIterator<VariantContext> it = vcfr.iterator();
        while(it.hasNext()) {
            VariantContext vc = it.next();
            String contig = vc.getContig(); // this should be the chromosome
            int pos = vc.getEnd(); // position of the variant in the genome
            String ref = vc.getReference().getDisplayString(); // get the reference Allele nucleotide
            String alt = vc.getAltAlleleWithHighestAlleleCount().getDisplayString();

            System.out.print("-------------------------------------------\n");
            System.out.print("chr" + contig + ":" + pos + ref + ">" + alt + "\t" + vc.getID() + "\n");

            /**/
            Map<String, Object> M = vc.getAttributes();
            for (String k: M.keySet()) {
                System.out.print(k + "\n");

            }
            /**/
            System.out.print("-------------------------------------------\n");
            break;
        }

        vcfr.close();
    }

}
