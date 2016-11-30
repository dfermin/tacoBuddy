package sampsonLab;

import org.apache.commons.cli.HelpFormatter;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.zip.GZIPInputStream;


/**
 * Created by dfermin on 11/30/16.
 */


public class globalFunctions {

    static public File inputVCF = null;
    static public File inputVCF_tabix = null;
    static File srcGFF3 = null;
    static String filterDOM = null;
    static String filterREC = null;
    static Set<String> genesDOM = null;
    static Set<String> genesREC = null;
    static double popFilter = 0.01; // default value
    static String tsOutputModel = null;
    static HashMap<String, Gene> REC_geneMap = null, DOM_geneMap = null;

    static public void parseCommandLineArgs(String[] args) throws IOException {
        File paramF = null;

        if( args[0].equalsIgnoreCase("-t") ) {
            writeTemplateInputFile();
            System.exit(0);
        }

        if( args[0].equalsIgnoreCase("-i") ) {
            paramF = new File(args[1]);
            if (!paramF.exists()) {
                System.err.print("\nERROR! Unable to find " + paramF.getCanonicalPath() + "\n");
                System.exit(0);
            }
            else {
                parseParamFile(paramF);
            }
        }
    }


    private static void parseParamFile(File paramF) throws IOException {
        FileReader fr = new FileReader(paramF);
        BufferedReader br = new BufferedReader(fr);

        String line;
        while((line = br.readLine()) != null) {

            if( line.contains("#") ) continue; // we skip lines that contain '#' character
            if( line.length() < 4 ) continue; // not enough data on this line to parse

            if(line.startsWith("inputVCF=")) {
                String vcf = line.substring(9);
                String tbi = vcf + ".tbi";
                inputVCF = new File(vcf);
                inputVCF_tabix = new File(tbi);
                continue;
            }

            if(line.startsWith("srcGFF3=")) {
                String gff3 = line.substring(8);
                srcGFF3 = new File(gff3);
                continue;
            }

            if(line.startsWith("genesDOM=")) {
                genesDOM = new HashSet<String>();
                String[] tmp = line.substring(9).split(",");
                for(String s: tmp) genesDOM.add(s.replaceAll("\\s*", ""));
            }
            if(line.startsWith("genesREC=")) {
                genesREC = new HashSet<String>();
                String[] tmp = line.substring(9).split(",");
                for(String s: tmp) genesREC.add(s.replaceAll("\\s*", ""));
            }

            if(line.startsWith("filterDOM=")) filterDOM = line.substring(10);
            if(line.startsWith("filterREC=")) filterREC = line.substring(10);

            if(line.startsWith("popFilter=")) popFilter = Double.parseDouble(line.substring(10));
            if(line.startsWith("tsOutputModel=")) tsOutputModel = line.substring(14);
        }
        br.close();

        errorCheck_paramFile_input();

        // If you go this far you had an acceptable input file. Report all copied values to STDERR
        System.err.print("\n-------------------------------------------------------------\n");
        if( !(null == inputVCF) )    System.err.print("inputVCF:     " + inputVCF.getCanonicalPath() + "\n");
        if( !(null == srcGFF3) )     System.err.print("source GFF:   " + srcGFF3.getCanonicalPath() + "\n");
        if( genesDOM.size() > 0 )  System.err.print("DOM genes:    " + genesDOM + "\n");
        if( genesREC.size() > 0 )  System.err.print("REC genes:    " + genesREC + "\n");
        if( filterDOM.length() > 0 ) System.err.print("DOM filters:  " + filterDOM + "\n");
        if( filterREC.length() > 0 ) System.err.print("REC filters:  " + filterREC + "\n");
        if( tsOutputModel.length() > 0 ) System.err.print("Output:       " + tsOutputModel + "\n");
        System.err.print("-------------------------------------------------------------\n\n");
    }



    private static void errorCheck_paramFile_input() {
        int score = 0;

        if( !inputVCF.exists() ) {
            System.err.println("\nERROR! Unable to find '" + inputVCF.getName() + "'\n");
            System.exit(0);
        }

        // This code assumes you have have a tabix *.tbi file to go with the input VCF file.
        if( !inputVCF_tabix.exists() ) {
            System.err.println(
                    "\nERROR! Unable to find tabix index file for '" + inputVCF_tabix.getName() + "'\n" +
                            "If you don't have it first create it with:\n" +
                            "\ttabix " + inputVCF.getName() + "\n"
            );
            System.exit(0);
        }

        if( !srcGFF3.exists() ) {
            System.err.println("\nERROR! Unable to find '" + srcGFF3.getName() + "'\n");
            System.exit(0);
        }

        score = 0;
        if( genesDOM.size() > 0 ) score++;
        if( genesREC.size() > 0 ) score++;
        if( score == 0 ) {
            System.err.println("\nERROR! You must provide a value for EITHER genesDOM or genesREC (or both) in the input file\n");
            System.exit(0);
        }

        score = 0;
        if( filterDOM.length() > 0 ) score++;
        if( filterREC.length() > 0 ) score++;
        if( score == 0 ) {
            System.err.println("\nERROR! You must provide a value for EITHER filterDOM or filterREC (or both) in the input file\n");
            System.exit(0);
        }
    }



    private static void writeTemplateInputFile() throws IOException {
        File templateF = new File("./tacoBuddy-inputTemplate.txt");
        FileWriter fw = new FileWriter(templateF);
        BufferedWriter bw = new BufferedWriter(fw);

        bw.write("# All lines containing a '#' are ignored.\n");

        bw.write("\n# This is the VCF file you want to analyze\ninputVCF=\n");
        bw.write("\n# This is the path to the gene coordinate file to use. Must be in GFF3 format.\nsrcGFF3=\n");
        bw.write("\n# List of DOMINANT genes to report results for\ngenesDOM=\n");
        bw.write("\n# List of RECESSIVE genes to report results for\ngenesREC=\n");
        bw.write("\n# Score filters to apply to the variants in DOMINANT genes\nfilterDOM=\n");
        bw.write("\n# Score filters to apply to the variants in RECESSIVE genes\nfilterREC=\n");
        bw.write("\npopFilter=0.01\n");
        bw.write("\n# The transcript model to report results for.\n#" +
                "Possible options are LT(longest transcript), MCT(most conserved transcript)\n" +
                "tsOutputModel=\n");
        bw.write("\n");
        bw.close();

        System.err.print("\nEdit template file: " + templateF.getCanonicalPath() + "\n");
    }



    public void parseGFF3() throws IOException {

        REC_geneMap = new HashMap<String, Gene>();
        DOM_geneMap = new HashMap<String, Gene>();
        BufferedReader br = null;

        System.err.print("\nParsing " + srcGFF3.getName() + "\n");

        if(srcGFF3.getName().endsWith(".gz")) {
            GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(srcGFF3.getAbsoluteFile()));
            br = new BufferedReader(new InputStreamReader(gzip));
        }
        else {
            FileReader fr = new FileReader(srcGFF3);
            br = new BufferedReader(fr);
        }

        String line;
        Gene curGene = null;
        while((line = br.readLine()) != null) {
            if(line.startsWith("#")) continue;

            curGene = new Gene(line);
            if(curGene.isValidGene()) {
                String k = curGene.getGeneName();
                if (genesDOM.contains(k)) DOM_geneMap.put(k, curGene);
                else if (genesREC.contains(k)) REC_geneMap.put(k, curGene);
            }
            curGene = null;
        }
        br.close();

        System.err.print(
                "DOM gene map size: " + DOM_geneMap.size() +
                "\nREC gene map size: " + REC_geneMap.size() + "\n");
    }
}