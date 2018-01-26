/**
 * Created by dfermin on 11/29/16.
 */

package sampsonLab;

import java.io.IOException;
import java.sql.SQLException;

import com.google.common.base.Joiner;


public class tacoBuddy {


    static public globalFunctions globals;
    static VCFParser the_vcf_parser;

    public static void main(String[] args) throws IOException, SQLException {

        // Prepare the globals object to handle user input
        globals = new globalFunctions();

        if( args.length < 1 ) {
            System.err.print("\nUSAGE: java -jar tacoBuddy.jar -i <input_file> or -t or -L <vcf.gz>\n" +
                    "\t-t to generate template input file\n" +
                    "\t-L list all the available INFO fields you can query or fetch from the given VCF.gz file\n\n");
            return;
        }

        globals.parseCommandLineArgs(args);
        globals.parseGFF3();
        globals.selectTranscriptModel();

        the_vcf_parser = new VCFParser(globals.inputVCF, globals.inputVCF_tabix);


        if(globals.featureSet.contains("EFF")) {
            double x = globals.required_min_sample_maf * 100;
            System.err.println("\nFilters based on EFF info fields:\n" +
                    "\t1. Synonymous SNPs excluded\n" +
                    "\t2. HIGH impact variants with Sample_MAF < " + x + "% will always be reported\n");
        }

        // Print header line
        String hdr = "#SAMPLE\tGenotype\tmodel\tTranscriptID\tGene_Name\tchrom\tpos\tdbsnpID\tREF\tALT\t" +
                     "ReadCounts (total)\tReadCounts (ref/alt)\tSAMPLE_MAF (as %)\t";
        hdr += Joiner.on("\t").join(globalFunctions.featureSet);
        System.out.println(hdr);

        // Record any variants that land within the exons of the genes stored in DOM and REC maps
        // If the genesDOM or genesREC sets are empty, we assume the user wants us to process all genes
        if(globals.queryMode.equalsIgnoreCase("exon")) {

            for(String geneType :  new String[] { "DOM", "REC"}) {

                if (geneType.equalsIgnoreCase("DOM")) {
                    if (globals.doAllGenes) the_vcf_parser.parseByExon(globals.ALL_geneMap, globals.filterDOM, geneType);
                    else if( null != globals.DOM_geneMap ) the_vcf_parser.parseByExon(globals.DOM_geneMap, globals.filterDOM, geneType);
                }

                if (geneType.equalsIgnoreCase("REC")) {
                    if (globals.doAllGenes) the_vcf_parser.parseByExon(globals.ALL_geneMap, globals.filterREC, geneType);
                    else if( null != globals.REC_geneMap ) the_vcf_parser.parseByExon(globals.REC_geneMap, globals.filterREC, geneType);
                }
            }
        }


        // Record any variants that land within the transcripts of the genes stored in DOM and REC maps
        if(globals.queryMode.equalsIgnoreCase("transcript")) {

            for(String geneType :  new String[] { "DOM", "REC"}) {

                if(geneType.equalsIgnoreCase("DOM")) {
                    if (globals.doAllGenes) the_vcf_parser.parseByTranscript(globals.ALL_geneMap, globals.filterDOM, geneType);
                    else if( null != globals.DOM_geneMap ) the_vcf_parser.parseByTranscript(globals.DOM_geneMap, globals.filterDOM, geneType);
                }

                if(geneType.equalsIgnoreCase("REC")) {
                    if (globals.doAllGenes) the_vcf_parser.parseByTranscript(globals.ALL_geneMap, globals.filterREC, geneType);
                    else if( null != globals.REC_geneMap ) the_vcf_parser.parseByTranscript(globals.REC_geneMap, globals.filterREC, geneType);
                }
            }
        }
    }

}
