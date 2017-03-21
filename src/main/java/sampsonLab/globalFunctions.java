package sampsonLab;

import com.google.common.base.Joiner;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.SetMultimap;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;


/**
 * Created by dfermin on 11/30/16.
 */


public class globalFunctions {

    static public File inputVCF = null;
    static public File inputVCF_tabix = null;
    static public File srcGFF3 = null;
    static public String filterDOM = null;
    static public String filterREC = null;
    static public Set<String> genesDOM = null;
    static public Set<String> genesREC = null;
    static public SortedSet<String> featureSet = null;
    static public SetMultimap<String, Transcript> REC_geneMap = null, DOM_geneMap = null;
    static public Map<String, String> allowedSitesMap = null;




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

            if( line.contains("GERP++") ) line = line.replaceAll("\\+\\+", "");

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
                String[] tmp = line.substring(9).split("[,\\t;\\s]");
                for(String s: tmp) genesDOM.add(s.replaceAll("\\s*", ""));
            }
            if(line.startsWith("genesREC=")) {
                genesREC = new HashSet<String>();
                String[] tmp = line.substring(9).split("[,\\t;\\s]");
                for(String s: tmp) genesREC.add(s.replaceAll("\\s*", ""));
            }

            if(line.startsWith("featureList=")) {
                featureSet = new TreeSet<String>();
                for(String s : line.substring(12).split("[,;\\s]+")) {
                    featureSet.add(s.toUpperCase()); // store all features as uppercase
                }
            }

            if(line.startsWith("allowedSites")) {
                if(null == allowedSitesMap) allowedSitesMap = new HashMap<String, String>();

                String tmp, tag;
                tmp = line.split("=")[0];
                if(tmp.endsWith("REC")) tag = "REC";
                else tag = "DOM";

                for(String s : line.substring(16).split("[,;\\s]+")) {
                    allowedSitesMap.put(s, tag);
                }
            }

            if(line.startsWith("filterDOM=")) filterDOM = reformat_filter(line.substring(10));
            if(line.startsWith("filterREC=")) filterREC = reformat_filter(line.substring(10));

        }
        br.close();

        errorCheck_paramFile_input();

        // If you go this far you had an acceptable input file. Report all copied values to STDERR
        System.err.print("\n-------------------------------------------------------------\n");
        if( !(null == inputVCF) ) System.err.print("inputVCF:     " + inputVCF.getCanonicalPath() + "\n");
        if( !(null == srcGFF3) ) System.err.print("source GFF:   " + srcGFF3.getCanonicalPath() + "\n");
        if( genesDOM.size() > 0 ) System.err.print("DOM genes:    " + genesDOM + "\n");
        if( genesREC.size() > 0 ) System.err.print("REC genes:    " + genesREC + "\n");
        if( filterDOM != null ) System.err.print("DOM filters:  " + filterDOM + "\n");
        if( filterREC != null ) System.err.print("REC filters:  " + filterREC + "\n");
        System.err.print("Selected output features: " + Joiner.on(", ").join(featureSet) + "\n");

        if( !(null == allowedSitesMap) && (allowedSitesMap.size() > 0) ) {
            System.err.print("Requested sites:\n");
            for(String k : allowedSitesMap.keySet()) { System.err.println("  " + allowedSitesMap.get(k) + "\t" + k); }
        }
        //if( tsOutputModel.length() > 0 ) System.err.print("Output:       " + tsOutputModel + "\n");
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
        if( (null != filterDOM) && (filterDOM.length() > 0) ) score++;
        if( (null != filterREC) && (filterREC.length() > 0) ) score++;
        if( score == 0 ) {
            System.err.println("\nERROR! You must provide a value for EITHER filterDOM or filterREC (or both) in the input file\n");
            System.exit(0);
        }
    }


    // Function reformats the user filter to make sure it will work with jexl
    // Specifically it converts: /[AD] =~ SIFT/ to ('A' =~ SIFT) or ('D' =~ SIFT)
    private static String reformat_filter(String origFilter) {
        String ret = "";
        if(!origFilter.contains("/")) ret = origFilter;
        else {
            ArrayList<String> AL = new ArrayList<String>();
            ArrayList<Integer> slashPos = new ArrayList<Integer>();
            for(int i = 0; i < origFilter.length(); i++) {
                if(origFilter.charAt(i) == '/') slashPos.add(i);
            }

            // slashPos must have an even number of entries since the '/' always comes as a pair
            for(int i = 0; i < slashPos.size() - 1; i += 2) {
                int j = i + 1;
                String tmp = origFilter.substring(slashPos.get(i)+1, slashPos.get(j));
                AL.add(tmp);
            }

            // AL now contains every instance of a perl-like or regex (example: ['ABCD'] =~ blah)
            // Reformat each one to be in compliance with the jexl api
            HashMap<String, String> jexlMap = new HashMap<String, String>();
            String tmp = "";
            for(String s : AL) {
                String key = s;
                s = s.trim().replaceAll("\\[", "").replaceAll("\\]", "").replaceAll("'", ""); // strip out []'
                int i = s.indexOf("=~");
                String comparison_suffix = s.substring(i).trim();
                String pattern = s.substring(0,i).trim();
                ArrayList<String> Jarray = new ArrayList<String>();
                for(String k : pattern.split("")) {
                    tmp = "('" + k + "' " + comparison_suffix + ")";
                    Jarray.add(tmp);
                }
                tmp = "(" + Joiner.on(" or ").join(Jarray) + ")";
                jexlMap.put(key, tmp); // k = original string, v = jexl compliant string
            }

            // Now replace the regex values in origFilter with the new values we've created in jexlMap
            ret = origFilter.replaceAll("/", "");
            for(String k : jexlMap.keySet()) {
                ret = ret.replace(k, jexlMap.get(k));
            }
        }
        return ret;
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
        bw.write("\n# List the variant information (ie: features) in the VCF file you want to report in the final output.\n" +
                        "# NOTE: This list _MUST_ contain all of the field names you use in 'filterDOM' and filterREC' below.\n" +
                        "# The entries here can be separated by tabs, spaces, commas or semicolons\n" +
                        "featureList=GENEINFO\n");
        bw.write("\n# Enter regex-like filters there for genes. For 'or' conditions surround the whole regex in forward slashes. Example: /['DT'] =~ SIFT/\n");
        bw.write("\n# Score filters to apply to the variants in DOMINANT genes\nfilterDOM=\n");
        bw.write("\n# Score filters to apply to the variants in RECESSIVE genes\nfilterREC=\n");
        bw.write("\n# Include data for these variants regardless of their filter scores\n" +
                "# Variant syntax: chromosome:position\n" +
                "# Multiple sites can be given as comma separated.\n#allowedSitesDOM=\n#allowedSitesREC=\n");
        bw.write("\n");
        bw.close();

        System.err.print("\nEdit template file: " + templateF.getCanonicalPath() + "\n");
    }



    public void parseGFF3() throws IOException {

        REC_geneMap = HashMultimap.create(); // k = geneName, v = list of transcripts
        DOM_geneMap = HashMultimap.create(); // k = geneName, v = list of transcripts
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

        Transcript curTranscript = null;
        int geneStart = -1, geneEnd = -1;
        char strand = '.';
        String geneID = null, geneType = null, geneName = null, chr = null, line = null;

        while((line = br.readLine()) != null) {
            if(line.startsWith("#")) continue;

            String[] data = line.split("\t");
            if(data[2].equalsIgnoreCase("gene")) {

                chr = data[0];
                geneStart = Integer.parseInt(data[3]);
                geneEnd = Integer.parseInt(data[4]);
                strand = data[6].charAt(0);


                String info[] = data[8].split(";");
                for (String s : info) {
                    if (s.startsWith("ID=")) geneID = s.substring(3);
                    if (s.startsWith("gene_type=")) geneType = s.substring(10);
                    if (s.startsWith("gene_name=")) geneName = s.substring(10);
                }
            }

            if(data[2].equalsIgnoreCase("transcript")) {

                // Check to see if curTranscript is null, if it isn't you need to store this variable
                // before you can continue;
                if(curTranscript != null) {
                    if(genesDOM.contains(geneName)) { DOM_geneMap.put(geneName, curTranscript); }
                    else if(genesREC.contains(geneName)) { REC_geneMap.put(geneName, curTranscript); }
                    curTranscript = null;
                }

                int start = Integer.parseInt(data[3]);
                int end = Integer.parseInt(data[4]);

                String localGeneName = null, transID = null;

                String info[] = data[8].split(";");
                for (String s : info) {
                    if (s.startsWith("ID=")) transID = s.substring(3);
                    if (s.startsWith("gene_name=")) localGeneName = s.substring(10);
                }

                if (!localGeneName.equalsIgnoreCase(geneName)) {
                    System.err.print("\nERROR! (transcript): " + localGeneName + " != " + geneName + " for current line of srcGFF3:\n" +
                            line + "\n\n");
                    System.exit(0);
                }

                // Create a new transcript entry
                curTranscript = new Transcript(chr,
                        geneStart, geneEnd,
                        strand, geneID, geneType, geneName,
                        start, end, transID);
            }


            if(data[2].equalsIgnoreCase("exon")) {
                int start = Integer.parseInt(data[3]);
                int end = Integer.parseInt(data[4]);

                String transID = null, exonID = null;
                int exonNum = 0;
                String info[] = data[8].split(";");
                for (String s : info) {
                    if (s.startsWith("Parent=")) transID = s.substring(7);
                    if (s.startsWith("exon_id=")) exonID = s.substring(8);
                    if (s.startsWith("exon_number=")) exonNum = Integer.parseInt(s.substring(12));
                }

                if(!transID.equalsIgnoreCase(curTranscript.getTranscriptID())) {
                    System.err.print("\nERROR! (exon): " + transID + " != " + curTranscript.getTranscriptID() + " for current line of srcGFF3:" +
                            line + "\n\n");
                    System.exit(0);
                }

                Exon E = new Exon(exonID,
                        curTranscript.getGeneName(),
                        curTranscript.getChrom(),
                        start, end, exonNum
                        );

                curTranscript.addExon(E);
                E = null;

            }
        }
        br.close();

        System.err.print(
                "DOM gene map size: " + DOM_geneMap.asMap().size() +
                "\nREC gene map size: " + REC_geneMap.asMap().size() + "\n");
    }
}
