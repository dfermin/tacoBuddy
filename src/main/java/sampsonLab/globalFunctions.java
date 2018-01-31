package sampsonLab;

import com.google.common.base.Joiner;

import com.google.common.collect.HashMultimap;

import com.google.common.collect.SetMultimap;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;


/**
 * Created by dfermin on 11/30/16.
 */


public class globalFunctions {

    static public File inputVCF = null;
    static public File inputVCF_tabix = null;
    static public File srcGFF3 = null;
    static public File TIMSfile = null;
    static public String vcfFormat = null; // either 'gotcloud' or 'gatk'
    static public String filterDOM = null;
    static public String filterREC = null;
    static public String queryMode = null;
    static public String outputTranscript = null;
    static public double required_min_sample_af = 0.1;
    static public boolean doAllGenes = false;
    static public Set<String> genesDOM = null;
    static public Set<String> genesREC = null;
    static public SortedSet<String> featureSet = null;
    static public SetMultimap<String, Transcript> REC_geneMap = null, DOM_geneMap = null, ALL_geneMap = null;
    static public TreeMap<String, Integer> allowedSites = null;



    /*----------------------------------------------------------------------------------------------------------------*/
    static public void parseCommandLineArgs(String[] args) throws IOException {
        File paramF = null;

        vcfFormat = "gotCloud";

        for(int j = 0; j < args.length; j++) {

            if( args[j].equalsIgnoreCase("-t") ) {
                writeTemplateInputFile();
                System.exit(0);
            }

            if( args[j].equalsIgnoreCase("-L") ) {
                String vcf = args[1];
                String tbi = vcf + ".tbi";
                inputVCF = new File(vcf);
                inputVCF_tabix = new File(tbi);
                VCFParser vcfp = new VCFParser(inputVCF, inputVCF_tabix);
                vcfp.printINFOfields();
                System.exit(0);
            }

            if( args[j].endsWith("gatk") ) {
                vcfFormat = "GATK";
            }

            if( args[j].equalsIgnoreCase("-i") ) {
                paramF = new File(args[(j+1)]);
                if (!paramF.exists()) {
                    System.err.print("\nERROR! Unable to find " + paramF.getCanonicalPath() + "\n");
                    System.exit(0);
                }
                else {
                    parseParamFile(paramF);
                }
            }
        }
    }


    /*----------------------------------------------------------------------------------------------------------------*/
    private static void parseParamFile(File paramF) throws IOException {
        FileReader fr = new FileReader(paramF);
        BufferedReader br = new BufferedReader(fr);

        genesDOM = new HashSet<String>();
        genesREC = new HashSet<String>();

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

            if(line.startsWith("min_sample_af=")) {
                required_min_sample_af = Double.valueOf(line.substring(15));
            }

            if(line.startsWith("queryMode=")) {
                queryMode = line.substring(10);
            }

            if(line.startsWith("genesDOM=")) {
               String[] tmp = line.substring(9).split("[,\\t;\\s]");
                for(String s: tmp) {
                    String s2 = s.replaceAll("\\s*", "");
                    if(s2.length() > 1) genesDOM.add(s2);
                }
            }
            if(line.startsWith("genesREC=")) {
                String[] tmp = line.substring(9).split("[,\\t;\\s]");
                for(String s: tmp) {
                    String s2 = s.replaceAll("\\s*", "");
                    if(s2.length() > 1) genesREC.add(s2);
                }
            }

            if(line.startsWith("transcriptModel=")) {
                outputTranscript=line.substring(16);
            }

            if(line.startsWith("TIMS=")) { // you only need this if the user elects to get the 'most conserved transcript'
                TIMSfile = new File(line.substring(5));
            }

            if(line.startsWith("featureList=")) {
                featureSet = new TreeSet<String>();
                for(String s : line.substring(12).split("[,;\\s]+")) {
                    featureSet.add(s.toUpperCase()); // store all features as uppercase
                }
            }

            // Sites to keep and report no matter what their filter scores are
            if(line.startsWith("allowedSites")) {
                if(null == allowedSites ) {
                    allowedSites = new TreeMap<String, Integer>();
                }
                for(String s : line.substring(13).split("[,;\\s]+")) {
                    if(s.length() > 5) allowedSites.put(s, 0);
                }
                if(allowedSites.isEmpty()) allowedSites = null;
            }

            if(line.startsWith("filterDOM=")) filterDOM = line.substring(10);
            if(line.startsWith("filterREC=")) filterREC = line.substring(10);

        }
        br.close();

        errorCheck_paramFile_input();


        // If you go this far you had an acceptable input file. Report all copied values to STDERR
        System.err.print("\n-------------------------------------------------------------\n");
        System.err.print("inputVCF:         " + inputVCF.getCanonicalPath() + "\n");
        System.err.print("source GFF:       " + srcGFF3.getCanonicalPath() + "\n");
        System.err.print("Query mode:       " + queryMode + "\n");
        System.err.print("Sample MAF:     < " + required_min_sample_af + "\n");

        System.err.print("Transcript model: " + outputTranscript + "\n");
        if( outputTranscript.equalsIgnoreCase("mostConserved") ) System.err.print("TIMS file:    " + TIMSfile.getName() + "\n");

        if( vcfFormat == "GATK" ) {
            System.err.println("Assuming GATK formatted VCF file");
        }

        if(null != ALL_geneMap) {
            System.err.print("Processing all genes in:  " + srcGFF3.getName() + "\n");
        }
        else {
            if (!genesDOM.isEmpty()) System.err.print("\nDOM genes:    " + genesDOM + "\n");
            if (!genesREC.isEmpty()) System.err.print("\nREC genes:    " + genesREC + "\n");
        }

        if (filterDOM != null) System.err.print("\nDOM filters:  " + filterDOM + "\n");
        if (filterREC != null) System.err.print("\nREC filters:  " + filterREC + "\n");

        System.err.print("\nSelected output features: " + Joiner.on(", ").join(featureSet) + "\n");

        if( !(null == allowedSites) && (allowedSites.size() > 0) ) {
            System.err.print("Requested sites:\n");
            for(String k : allowedSites.keySet()) {
                System.err.println("\t" + k);
            }
        }
        System.err.print("-------------------------------------------------------------\n\n");
    }


    /*----------------------------------------------------------------------------------------------------------------*/
    private static void errorCheck_paramFile_input() {
        int score = 0;

        if( !inputVCF.exists() ) {
            System.err.println("\nERROR! Unable to find '" + inputVCF.getName() + "'\n");
            System.exit(1);
        }

        // This code assumes you have have a tabix *.tbi file to go with the input VCF file.
        if( !inputVCF_tabix.exists() ) {
            System.err.println(
                    "\nERROR! Unable to find tabix index file for '" + inputVCF_tabix.getName() + "'\n" +
                            "If you don't have it first create it with:\n" +
                            "\ttabix " + inputVCF.getName() + "\n"
            );
            System.exit(1);
        }

        if( !srcGFF3.exists() ) {
            System.err.println("\nERROR! Unable to find '" + srcGFF3.getName() + "'\n");
            System.exit(1);
        }

        if( outputTranscript.equalsIgnoreCase("mostConserved") &&
                ( (null == TIMSfile) || !TIMSfile.exists() ) ) {
            System.err.println("\nERROR! You asked for the 'mostConserved' transript model output but I can't find the TIMS file\n");
            System.exit(1);
        }

        if( null == queryMode ) {
            System.err.println("\nERROR! queryMode is not set. Must be either 'exon' or 'transcript'\n");
            System.exit(1);
        }

        if( !queryMode.equalsIgnoreCase("transcript") &&
            !queryMode.equalsIgnoreCase("exon")) {
            System.err.println("\nERROR! queryMode must be either 'exon' or 'transcript'\n");
            System.exit(1);
        }


        // If the user has not specified any DOM or REC genes then initialize ALL_geneMap
        score = 0;
        if( (null != genesDOM) && (genesDOM.size() > 0) ) score++;
        if( (null != genesREC) && (genesREC.size() > 0) ) score++;
        if(score == 0) {
            ALL_geneMap = HashMultimap.create();
            doAllGenes = true;
        }


//        score = 0;
//        if( (null != filterDOM) && (filterDOM.length() > 0) && (genesDOM.size() > 0) ) score++;
//        if( (null != filterREC) && (filterREC.length() > 0)  && (genesREC.size() > 0) ) score++;
//        if( score == 0 ) {
//            System.err.println("\nERROR! You must provide a value for EITHER filterDOM or filterREC (or both) in the input file\n");
//            System.exit(1);
//        }

        if( (genesDOM.size() > 0) && (filterDOM.length() == 0) ) {
            System.err.println("\nERROR! You provided " + genesDOM.size() + " DOMINANT model genes but no filter string for them.\n");
            System.exit(1);
        }

        if( (genesREC.size() > 0) && (filterREC.length() == 0) ) {
            System.err.println("\nERROR! You provided " + genesREC.size() + " RECESSIVE model genes but no filter string for them.\n");
            System.exit(1);
        }

        checkFilterStrings(filterDOM);
        checkFilterStrings(filterREC);
    }

    /*----------------------------------------------------------------------------------------------------------------*/
    // Function parses the given filter string to make sure all of the features in the string are also
    // listed in the featuresList object. Any missing features will be added.
    public static void checkFilterStrings(String filter_str) {

        ArrayList<String> bits = new ArrayList<String>();

        if(null == filter_str) return;

        // Strip the string of parentheses
        filter_str = filter_str.replaceAll("[\\(\\)'\\/]+", "");

        // First split the string up by 'and' & 'or' statements
        String[] parts = filter_str.split("and|or");

        for(String p : parts) {

            // Now split this string up by regex symbols
            String[] parts2 = p.split("[=~<>\\s]+");
            for(String p2 : parts2) {

                if(Pattern.matches("^[\\d\\.-]+$", p2)) continue; // just a number
                if(p2.contains("[") || p2.contains("]") ) continue; // a regex box
                if(p2.length() < 3) continue;
                if(p2.equalsIgnoreCase("sample_maf")) continue; // this is added automatically
                bits.add(p2);
            }
        }

        // 'bits' now contains only words/terms to filter on.
        // Curate the data in 'bits' removing any special values that are not going to match
        // a VCF INFO field but are instead computed internally by this program.
        ArrayList<String> cleanBits = new ArrayList<String>();
        for(String b : bits) {
            if(b.equalsIgnoreCase("HIGH")) continue;
            if(b.equalsIgnoreCase("MODERATE")) continue;
            if(b.equalsIgnoreCase("LOW")) continue;

            cleanBits.add(b);
        }

        // Go through cleanBits and make sure they are in the final output feature list
        for(String b : cleanBits) {
            if( !featureSet.contains(b.toUpperCase()) ) featureSet.add(b);
        }
    }

    /*----------------------------------------------------------------------------------------------------------------*/
    private static void writeTemplateInputFile() throws IOException {
        File templateF = new File("./tacoBuddy-inputTemplate.txt");
        FileWriter fw = new FileWriter(templateF);
        BufferedWriter bw = new BufferedWriter(fw);

        bw.write("# All lines containing a '#' are ignored.\n");

        bw.write("\n# This is the VCF file you want to analyze\ninputVCF=\n");
        bw.write("\n# This is the path to the gene coordinate file to use. Must be in GFF3 format.\nsrcGFF3=\n");
        bw.write("\n# List of DOMINANT genes to report results for\ngenesDOM=\n");
        bw.write("\n# List of RECESSIVE genes to report results for\ngenesREC=\n");
        bw.write("\n# Specify the minimum Sample Allele Frequecy (AF) that a variant call must have in order to be reported");
        bw.write("\nmin_sample_af=0.05");
        bw.write("\n# Specify query method for reporting variant calls." +
                      "\n# 'transcript' = look for variants within the boundaries of a transcript" +
                      "\n# 'exon' = look for variants within the coding boundaries of exons" +
                      "\nqueryMode=transcript\n");

        bw.write("\n# List the variant information (ie: features) in the VCF file you want to report in the final output." +
                      "\n# NOTE: This list _MUST_ contain all of the field names you use in 'filterDOM' and filterREC' below." +
                      "\n# The entries here can be separated by tabs, spaces, commas or semicolons" +
                      "\nfeatureList=EFF, IS_LOF, EXAC_KG_AF_POPMAX, ONEKG_AF\n");

        bw.write("\n# Enter regex-like filters there for genes. " +
                      "\n# The syntax is **CRITICAL** here!" +
                      "\n# You must write your text-based expressions as: 'vcf_field_name' =~ 'query_term' " +
                      "\n# Numerical filters are written as 'filter_field' <= 'numeric_cutoff' " +
                      "\n# Some default filters are given below as examples of proper syntax usage." +
                      "\n# These examples assume that dbNSFP, EFF, and ESP were used to annotate the VCF\n");
        bw.write("\n# Score filters to apply to the variants in DOMINANT genes" +
                      "\nfilterDOM=SAMPLE_AF < 0.1 and (ESP_MAX_AA_EA < 0.005 and (( (POLYPHEN2_HVAR =~ /D/) + (MUTATIONTASTER =~ /[AD]/) + (SIFT =~ /D/) ) >= 2)) or (ESP_MAX_AA_EA < 0.005 and IS_LOF)\n");

        bw.write("\n# Score filters to apply to the variants in RECESSIVE genes" +
                      "\nfilterREC=SAMPLE_AF < 0.1 and (ESP_MAX_AA_EA < 0.01 and (( (POLYPHEN2_HVAR =~ /D/) + (MUTATIONTASTER =~ /[AD]/) + (SIFT =~ /D/) ) >= 2)) or (ESP_MAX_AA_EA < 0.01 and IS_LOF)\n");

        bw.write("\n# Include data for these variants regardless of their filter scores\n" +
                "# Variant syntax: chromosome:position\n" +
                "# Multiple sites can be given as comma separated." +
                "\n#allowedSites=\n");

        bw.write("\n# Specify which transcript model you want as output. The options are: all, longest, or mostConserved" +
                "\n# mostConserved is defined by the TIMs score available at this site: http://glom.sph.umich.edu/GIMS" +
                "\n# If you go with the mostConserved option, you *MUST* also provide the TIMS score file from the above website." +
                "\n# The 'all' option will return ALL transcripts associated with the gene. This option is more for debugging." +
                "\ntranscriptModel=longest\n");
        bw.write("\n# You only need to provide this file if you chose 'mostConserved' for the 'transcriptModel' option\n# TIMS=\n");
        bw.write("\n");
        bw.close();

        System.err.print("\nEdit template file: " + templateF.getCanonicalPath() + "\n");
    }


    /*----------------------------------------------------------------------------------------------------------------*/
    public void parseGFF3() throws IOException {

        if( !doAllGenes ) {
            if (this.genesREC.size() > 0) REC_geneMap = HashMultimap.create(); // k = geneName, v = list of transcripts
            if (this.genesDOM.size() > 0) DOM_geneMap = HashMultimap.create(); // k = geneName, v = list of transcripts
        }


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
                    if (s.startsWith("ID=")) geneID = s.substring(3).replaceAll("\\..+$", "");  // remove version number from ENSG id
                    if (s.startsWith("gene_type=")) geneType = s.substring(10);
                    if (s.startsWith("gene_name=")) geneName = s.substring(10);
                }
            }

            if(data[2].equalsIgnoreCase("transcript")) {

                // Check to see if curTranscript is null, if it isn't you need to store this variable first
                // before you can continue;
                if(curTranscript != null) {
                    curTranscript.calcCDSlength();
                    if(genesDOM.contains(curTranscript.getGeneName()))
                        DOM_geneMap.put(curTranscript.getGeneName(), curTranscript);

                    else if(genesREC.contains(curTranscript.getGeneName()))
                        REC_geneMap.put(curTranscript.getGeneName(), curTranscript);

                    else if( null != ALL_geneMap )
                        ALL_geneMap.put(curTranscript.getGeneName(), curTranscript);
                }
                curTranscript = null;

                int start = Integer.parseInt(data[3]);
                int end = Integer.parseInt(data[4]);

                String localGeneName = null, transID = null;

                String info[] = data[8].split(";");
                for (String s : info) {
                    if (s.startsWith("ID=")) transID = s.substring(3).replaceAll("\\..+$", ""); // remove version number from ENST id
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
                        strand, geneID, geneType, localGeneName,
                        start, end, transID);
            }


            if(data[2].equalsIgnoreCase("exon")) {
                int start = Integer.parseInt(data[3]);
                int end = Integer.parseInt(data[4]);

                String transID = null, exonID = null;
                int exonNum = 0;
                String info[] = data[8].split(";");
                for (String s : info) {
                    if (s.startsWith("Parent=")) transID = s.substring(7).replaceAll("\\..+$", "");
                    if (s.startsWith("exon_id=")) exonID = s.substring(8).replaceAll("\\..+$", "");
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

        if(doAllGenes) System.err.println("Gene map size: " + ALL_geneMap.asMap().size());
        else {
            if (this.genesDOM.size() > 0) System.err.println("DOM gene map size: " + DOM_geneMap.asMap().size());
            if (this.genesREC.size() > 0) System.err.println("REC gene map size: " + REC_geneMap.asMap().size());
        }
    }


    /*-----------------------------------------------------------------------------------------------------------------
    ** Function chooses which transcripts to report in the final output based upon the outputTranscript variable
     */
    public void selectTranscriptModel() throws IOException {

        HashMap<String, Transcript> gene2transMap = null;
        Transcript best_TS = null;

        if(this.outputTranscript.equalsIgnoreCase("all")) return; // don't filter transcript

        if(this.outputTranscript.equalsIgnoreCase("mostConserved")) {
            this.parseTIMS();

            double lowest_tims;
            gene2transMap = new HashMap<String, Transcript>();

            if(doAllGenes) {
                // iterate over all the genes in ALL_geneMap
                // Keep the transcript with the *LOWEST* tims score
                for (String gene : ALL_geneMap.keys()) {
                    lowest_tims = 10000;
                    best_TS = null;

                    for (Transcript ts : ALL_geneMap.get(gene)) {
                        if (ts.getTims() < lowest_tims) {
                            lowest_tims = ts.getTims();
                            best_TS = ts;
                        }
                    }

                    if (lowest_tims < 10000) { // we found at least 1 transcript with a low TIMS score
                        gene2transMap.put(gene, best_TS);
                    } else { // all of the transcripts have the same TIMS score so keep the longest transcript
                        int tsLen = 0;
                        best_TS = null;
                        for (Transcript k : ALL_geneMap.get(gene)) {
                            if (k.getTsLen() > tsLen) {
                                tsLen = k.getTsLen();
                                best_TS = k;
                            }
                        }
                        gene2transMap.put(gene, best_TS);
                    }
                }
                ALL_geneMap.clear();
                for (String k : gene2transMap.keySet()) ALL_geneMap.put(k, gene2transMap.get(k));
                gene2transMap.clear();
            }
            else {
                // Iterate over all the genes in DOM_geneMap
                // Keep the transcript with the *LOWEST* tims score
                for (String gene : DOM_geneMap.keys()) {
                    lowest_tims = 10000;
                    best_TS = null;

                    for (Transcript ts : DOM_geneMap.get(gene)) {
                        if (ts.getTims() < lowest_tims) {
                            lowest_tims = ts.getTims();
                            best_TS = ts;
                        }
                    }

                    if (lowest_tims < 10000) { // we found at least 1 transcript with a low TIMS score
                        gene2transMap.put(gene, best_TS);
                    } else { // all of the transcripts have the same TIMS score so keep the longest transcript
                        int tsLen = 0;
                        best_TS = null;
                        for (Transcript k : DOM_geneMap.get(gene)) {
                            if (k.getTsLen() > tsLen) {
                                tsLen = k.getTsLen();
                                best_TS = k;
                            }
                        }
                        gene2transMap.put(gene, best_TS);
                    }
                }
                DOM_geneMap.clear();
                for (String k : gene2transMap.keySet()) DOM_geneMap.put(k, gene2transMap.get(k));
                gene2transMap.clear();


                // Iterate over all the genes in REC_geneMap
                // Keep the transcript with the *LOWEST* tims score
                for (String gene : REC_geneMap.keys()) {
                    lowest_tims = 10000;
                    best_TS = null;

                    for (Transcript ts : REC_geneMap.get(gene)) {
                        if (ts.getTims() < lowest_tims) {
                            lowest_tims = ts.getTims();
                            best_TS = ts;
                        }
                    }

                    if (lowest_tims < 10000) {
                        gene2transMap.put(gene, best_TS);
                    } else { // all of the transcripts have the same TIMS score so keep the longest transcript
                        int tsLen = 0;
                        best_TS = null;
                        for (Transcript k : REC_geneMap.get(gene)) {
                            if (k.getTsLen() > tsLen) {
                                tsLen = k.getTsLen();
                                best_TS = k;
                            }
                        }
                        gene2transMap.put(gene, best_TS);
                    }
                }
                REC_geneMap.clear();
                for (String k : gene2transMap.keySet()) REC_geneMap.put(k, gene2transMap.get(k));
                gene2transMap.clear();

            } // end else

        } //-------------- End if over outputTranscript == mostConserved ------------------------




        if(outputTranscript.equalsIgnoreCase("longest")) {
            int tsLen = 0;
            best_TS = null;

            gene2transMap = new HashMap<String, Transcript>();

            if(doAllGenes) {
                for (String gene : ALL_geneMap.keys()) {
                    tsLen = 0;
                    best_TS = null;
                    for (Transcript ts : ALL_geneMap.get(gene)) {
                        if (ts.getTsLen() > tsLen) {
                            tsLen = ts.getTsLen();
                            best_TS = ts;
                        }
                    }
                    gene2transMap.put(gene, best_TS);
                }
                ALL_geneMap.clear();
                for (String k : gene2transMap.keySet()) ALL_geneMap.put(k, gene2transMap.get(k));
                gene2transMap.clear();
            }
            else {
                if (null != DOM_geneMap) {
                    for (String gene : DOM_geneMap.keys()) {
                        tsLen = 0;
                        best_TS = null;
                        for (Transcript ts : DOM_geneMap.get(gene)) {
                            if (ts.getTsLen() > tsLen) {
                                tsLen = ts.getTsLen();
                                best_TS = ts;
                            }
                        }
                        gene2transMap.put(gene, best_TS);
                    }
                    DOM_geneMap.clear();
                    for (String k : gene2transMap.keySet()) DOM_geneMap.put(k, gene2transMap.get(k));
                    gene2transMap.clear();
                }


                if (null != REC_geneMap) {
                    for (String gene : REC_geneMap.keys()) {
                        tsLen = 0;
                        best_TS = null;
                        for (Transcript ts : REC_geneMap.get(gene)) {
                            if (ts.getTsLen() > tsLen) {
                                tsLen = ts.getTsLen();
                                best_TS = ts;
                            }
                        }
                        gene2transMap.put(gene, best_TS);
                    }
                    REC_geneMap.clear();
                    for (String k : gene2transMap.keySet()) REC_geneMap.put(k, gene2transMap.get(k));
                    gene2transMap.clear();
                }
            } // end else
        } // end longest transcript
    }

    /*------------------------------------------------------------------------------------------------------------------
    ** Function parses the TIM file aquried from http://glom.sph.umich.edu/GIMS/
     */
    public void parseTIMS() throws IOException {

        System.err.println("Parsing " + this.TIMSfile.getName());

        FileReader fr = new FileReader(this.TIMSfile);
        BufferedReader br = new BufferedReader(fr);
        String line = null, geneID = null, transID = null;
        double tims = -1;

        while((line = br.readLine()) != null) {
            if(line.startsWith("GENCODE")) continue;

            String lineData[] = line.split("\\s+");
            geneID = lineData[1]; // gene Symbol
            transID = lineData[2].replaceAll("\\..+$", "");
            tims = Double.valueOf( lineData[3] );

            if( (null != DOM_geneMap) && (DOM_geneMap.containsKey(geneID)) ) {
                for(Transcript ts : DOM_geneMap.get(geneID)) {
                    if(ts.getTranscriptID().equalsIgnoreCase(transID)) {
                        ts.setTIMS(tims);
                        DOM_geneMap.put(geneID, ts);
                    }
                }
            }


            if( (null != REC_geneMap) && (REC_geneMap.containsKey(geneID)) ) {
                for(Transcript ts : REC_geneMap.get(geneID)) {
                    if(ts.getTranscriptID().equalsIgnoreCase(transID)) {
                        ts.setTIMS(tims);
                        REC_geneMap.put(geneID, ts);
                    }
                }
            }

        }
        br.close();
    }


    /*------------------------------------------------------------------------------------------------------------------
    ** Function takes in an Object that is an 'ArrayList' object and
    ** 1) makes it non-redundant
    ** 2) removes any NULLs *IF* their are other values in the arrayList
     */
    public static String  arrayList2String(Object o, String nullStr) {

        // Assign the contents of the arrayList to a hashSet (thus making them unique)
        Set<String> hs = new HashSet<String>();
        hs.addAll((ArrayList<String>) o);

        // Store the hashSet back into an arrayList
        ArrayList<String> al = new ArrayList<String>();
        al.addAll(hs);

        // If the arrayList as more than one value and one of them is 'nullStr', remove 'nullStr'
        if ((al.size() > 1) && (al.contains(nullStr))) {
            al.remove(nullStr);
        }

        String ret = String.join(",", al);
        return(ret);
    }


}
