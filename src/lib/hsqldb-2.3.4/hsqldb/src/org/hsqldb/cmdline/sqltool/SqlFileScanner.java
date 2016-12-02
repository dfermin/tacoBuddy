/* The following code was generated by JFlex 1.4.3 on 5/21/15 9:30 PM */

/* Copyright (c) 2001-2015, The HSQL Development Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * Neither the name of the HSQL Development Group nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL HSQL DEVELOPMENT GROUP, HSQLDB.ORG,
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


/* @(#)$Id: SqlFileScanner.java 5472 2015-05-22 23:00:58Z fredt $ */

package org.hsqldb.cmdline.sqltool;

import java.io.PrintStream;
import org.hsqldb.lib.FrameworkLogger;


/**
 * This class is a scanner generated by 
 * <a href="http://www.jflex.de/">JFlex</a> 1.4.3
 * on 5/21/15 9:30 PM from the specification file
 * <tt>/home/blaine/hsqldb/src/org/hsqldb/cmdline/sqltool.flex</tt>
 */
public class SqlFileScanner implements TokenSource {

  /** This character denotes the end of file */
  public static final int YYEOF = -1;

  /** initial size of the lookahead buffer */
  private static final int ZZ_BUFFERSIZE = 16384;

  /** lexical states */
  public static final int SPECIAL = 12;
  public static final int SQL_DOUBLE_QUOTED = 8;
  public static final int SQL_SINGLE_QUOTED = 6;
  public static final int GOBBLE = 10;
  public static final int RAW = 4;
  public static final int SQL = 2;
  public static final int YYINITIAL = 0;
  public static final int EDIT = 16;
  public static final int PL = 14;
  public static final int PROMPT_CHANGE_STATE = 20;
  public static final int MACRO = 18;

  /**
   * ZZ_LEXSTATE[l] is the state in the DFA for the lexical state l
   * ZZ_LEXSTATE[l+1] is the state in the DFA for the lexical state l
   *                  at the beginning of a line
   * l is of the form l = 2*k, k a non negative integer
   */
  private static final int ZZ_LEXSTATE[] = { 
     0,  0,  1,  2,  3,  3,  4,  4,  5,  5,  6,  6,  7,  7,  8,  8, 
     9,  9, 10, 10, 11, 11
  };

  /** 
   * Translates characters to character classes
   */
  private static final String ZZ_CMAP_PACKED = 
    "\11\0\1\6\1\2\1\0\1\6\1\1\22\0\1\6\1\0\1\3"+
    "\4\0\1\34\2\0\1\5\2\0\1\7\1\33\1\4\12\0\1\32"+
    "\1\10\5\0\1\20\1\11\1\16\1\26\1\12\1\22\1\13\1\0"+
    "\1\14\1\0\1\27\1\30\1\0\1\15\1\24\1\25\1\0\1\17"+
    "\1\0\1\21\1\23\6\0\1\31\4\0\1\20\1\11\1\16\1\26"+
    "\1\12\1\22\1\13\1\0\1\14\1\0\1\27\1\30\1\0\1\15"+
    "\1\24\1\25\1\0\1\17\1\0\1\21\1\23\uff8a\0";

  /** 
   * Translates characters to character classes
   */
  private static final char [] ZZ_CMAP = zzUnpackCMap(ZZ_CMAP_PACKED);

  /** 
   * Translates DFA states to action switch labels.
   */
  private static final int [] ZZ_ACTION = zzUnpackAction();

  private static final String ZZ_ACTION_PACKED_0 =
    "\14\0\1\1\2\2\1\3\1\4\1\5\1\6\1\7"+
    "\1\10\3\1\1\11\1\12\1\13\2\14\1\15\2\13"+
    "\1\16\1\17\2\20\1\13\1\0\2\21\1\0\1\22"+
    "\1\23\1\22\1\24\2\25\1\22\2\26\2\22\2\27"+
    "\2\30\2\31\2\32\7\0\2\33\10\0\2\34\2\35"+
    "\1\0\2\36\1\37\3\0\1\40\1\41\3\0\2\42"+
    "\42\0";

  private static int [] zzUnpackAction() {
    int [] result = new int[129];
    int offset = 0;
    offset = zzUnpackAction(ZZ_ACTION_PACKED_0, offset, result);
    return result;
  }

  private static int zzUnpackAction(String packed, int offset, int [] result) {
    int i = 0;       /* index in packed string  */
    int j = offset;  /* index in unpacked array */
    int l = packed.length();
    while (i < l) {
      int count = packed.charAt(i++);
      int value = packed.charAt(i++);
      do result[j++] = value; while (--count > 0);
    }
    return j;
  }


  /** 
   * Translates a state to a row index in the transition table
   */
  private static final int [] ZZ_ROWMAP = zzUnpackRowMap();

  private static final String ZZ_ROWMAP_PACKED_0 =
    "\0\0\0\35\0\72\0\127\0\164\0\221\0\256\0\313"+
    "\0\350\0\u0105\0\u0122\0\u013f\0\u015c\0\u0179\0\u015c\0\u015c"+
    "\0\u0196\0\u015c\0\u01b3\0\u01d0\0\u015c\0\u01ed\0\u020a\0\u0227"+
    "\0\u015c\0\u015c\0\u015c\0\u0244\0\u015c\0\u015c\0\u0261\0\u027e"+
    "\0\u015c\0\u015c\0\u029b\0\u015c\0\u02b8\0\u02d5\0\u02f2\0\u015c"+
    "\0\u030f\0\u032c\0\u0349\0\u0366\0\u015c\0\u0383\0\u015c\0\u015c"+
    "\0\u03a0\0\u015c\0\u03bd\0\u03da\0\u03f7\0\u015c\0\u0414\0\u015c"+
    "\0\u0431\0\u015c\0\u044e\0\u015c\0\u046b\0\u0488\0\u04a5\0\u04c2"+
    "\0\u04df\0\u04fc\0\u02b8\0\u0519\0\u015c\0\u0536\0\u0553\0\u0570"+
    "\0\u058d\0\u05aa\0\u05c7\0\u05e4\0\u0601\0\u061e\0\u015c\0\u063b"+
    "\0\u015c\0\u0658\0\u0675\0\u015c\0\u015c\0\u0692\0\u06af\0\u06cc"+
    "\0\u015c\0\u015c\0\u06e9\0\u0706\0\u0723\0\u0740\0\u015c\0\u075d"+
    "\0\u077a\0\u0797\0\u07b4\0\u07d1\0\u07ee\0\u080b\0\u0828\0\u0845"+
    "\0\u0862\0\u087f\0\u089c\0\u08b9\0\u08d6\0\u08f3\0\u0910\0\u092d"+
    "\0\u094a\0\u0967\0\u0984\0\u09a1\0\u09be\0\u09db\0\u09f8\0\u0a15"+
    "\0\u0a32\0\u0a4f\0\u0a6c\0\u0a89\0\u0aa6\0\u0ac3\0\u0ae0\0\u0afd"+
    "\0\u0b1a";

  private static int [] zzUnpackRowMap() {
    int [] result = new int[129];
    int offset = 0;
    offset = zzUnpackRowMap(ZZ_ROWMAP_PACKED_0, offset, result);
    return result;
  }

  private static int zzUnpackRowMap(String packed, int offset, int [] result) {
    int i = 0;  /* index in packed string  */
    int j = offset;  /* index in unpacked array */
    int l = packed.length();
    while (i < l) {
      int high = packed.charAt(i++) << 16;
      result[j++] = high | packed.charAt(i++);
    }
    return j;
  }

  /** 
   * The transition table of the DFA
   */
  private static final int [] ZZ_TRANS = zzUnpackTrans();

  private static final String ZZ_TRANS_PACKED_0 =
    "\1\15\1\16\1\17\1\20\1\21\1\22\1\23\1\24"+
    "\1\25\1\26\4\15\1\27\7\15\1\30\2\15\1\31"+
    "\1\32\1\15\1\20\1\33\1\34\1\35\1\36\1\37"+
    "\2\33\1\40\1\41\23\33\1\42\1\33\1\43\1\44"+
    "\1\36\1\37\1\33\1\45\1\40\1\41\23\33\1\42"+
    "\1\46\1\47\1\50\3\46\1\4\24\46\1\51\1\46"+
    "\34\52\1\53\3\54\1\55\31\54\1\7\1\56\1\57"+
    "\32\7\1\60\1\61\1\62\1\60\1\63\2\60\1\64"+
    "\26\60\1\65\1\66\1\60\1\63\2\60\1\64\26\60"+
    "\1\67\1\70\33\60\1\71\1\72\1\60\1\63\2\60"+
    "\1\64\25\60\1\0\1\73\1\74\71\0\1\17\37\0"+
    "\1\75\35\0\1\23\26\0\1\24\2\0\32\24\12\0"+
    "\1\76\41\0\1\77\27\0\1\100\24\0\1\35\37\0"+
    "\1\101\36\0\1\102\27\0\1\44\33\0\1\43\1\44"+
    "\3\0\1\103\26\0\1\46\1\47\1\50\32\46\2\0"+
    "\1\50\32\0\1\46\1\104\1\105\3\46\1\51\1\46"+
    "\1\106\24\46\34\52\35\0\1\60\3\54\1\0\31\54"+
    "\2\0\1\57\34\0\1\62\37\0\1\107\36\0\1\110"+
    "\27\0\1\66\34\0\1\70\34\0\1\72\34\0\1\74"+
    "\32\0\5\75\1\111\27\75\13\0\1\112\33\0\1\113"+
    "\40\0\1\114\16\0\5\101\1\115\27\101\1\102\1\116"+
    "\1\117\32\102\2\0\1\105\32\0\1\46\1\120\1\121"+
    "\3\46\1\106\26\46\5\107\1\122\27\107\1\110\1\123"+
    "\1\124\32\110\4\75\1\125\1\111\27\75\14\0\1\126"+
    "\40\0\1\127\44\0\1\130\4\0\4\101\1\131\1\115"+
    "\27\101\2\0\1\117\34\0\1\121\32\0\4\107\1\132"+
    "\1\122\27\107\2\0\1\124\47\0\1\133\40\0\1\134"+
    "\33\0\1\135\15\0\1\136\1\137\3\0\1\133\40\0"+
    "\1\140\41\0\1\141\17\0\1\137\40\0\1\142\40\0"+
    "\1\133\30\0\1\142\13\0\1\143\1\0\1\144\1\145"+
    "\32\0\1\146\30\0\1\147\34\0\1\150\1\151\31\0"+
    "\1\152\25\0\1\153\52\0\1\154\26\0\1\155\34\0"+
    "\1\156\24\0\1\153\10\0\1\157\33\0\1\160\45\0"+
    "\1\161\26\0\1\162\25\0\1\163\34\0\1\164\42\0"+
    "\1\165\30\0\1\166\45\0\1\167\35\0\1\170\21\0"+
    "\1\171\45\0\1\172\40\0\1\173\27\0\1\174\23\0"+
    "\1\175\37\0\1\175\37\0\1\176\33\0\1\171\15\0"+
    "\1\175\1\136\1\137\32\175\16\0\1\177\30\0\1\200"+
    "\30\0\1\201\34\0\1\201\13\0\1\143\2\0\1\145"+
    "\7\0";

  private static int [] zzUnpackTrans() {
    int [] result = new int[2871];
    int offset = 0;
    offset = zzUnpackTrans(ZZ_TRANS_PACKED_0, offset, result);
    return result;
  }

  private static int zzUnpackTrans(String packed, int offset, int [] result) {
    int i = 0;       /* index in packed string  */
    int j = offset;  /* index in unpacked array */
    int l = packed.length();
    while (i < l) {
      int count = packed.charAt(i++);
      int value = packed.charAt(i++);
      value--;
      do result[j++] = value; while (--count > 0);
    }
    return j;
  }


  /* error codes */
  private static final int ZZ_UNKNOWN_ERROR = 0;
  private static final int ZZ_NO_MATCH = 1;
  private static final int ZZ_PUSHBACK_2BIG = 2;

  /* error messages for the codes above */
  private static final String ZZ_ERROR_MSG[] = {
    "Unkown internal scanner error",
    "Error: could not match input",
    "Error: pushback value was too large"
  };

  /**
   * ZZ_ATTRIBUTE[aState] contains the attributes of state <code>aState</code>
   */
  private static final int [] ZZ_ATTRIBUTE = zzUnpackAttribute();

  private static final String ZZ_ATTRIBUTE_PACKED_0 =
    "\14\0\1\11\1\1\2\11\1\1\1\11\2\1\1\11"+
    "\3\1\3\11\1\1\2\11\2\1\2\11\1\1\1\11"+
    "\1\1\1\0\1\1\1\11\1\0\3\1\1\11\1\1"+
    "\2\11\1\1\1\11\3\1\1\11\1\1\1\11\1\1"+
    "\1\11\1\1\1\11\7\0\1\1\1\11\10\0\1\1"+
    "\1\11\1\1\1\11\1\0\1\1\2\11\3\0\2\11"+
    "\3\0\1\1\1\11\42\0";

  private static int [] zzUnpackAttribute() {
    int [] result = new int[129];
    int offset = 0;
    offset = zzUnpackAttribute(ZZ_ATTRIBUTE_PACKED_0, offset, result);
    return result;
  }

  private static int zzUnpackAttribute(String packed, int offset, int [] result) {
    int i = 0;       /* index in packed string  */
    int j = offset;  /* index in unpacked array */
    int l = packed.length();
    while (i < l) {
      int count = packed.charAt(i++);
      int value = packed.charAt(i++);
      do result[j++] = value; while (--count > 0);
    }
    return j;
  }

  /** the input device */
  private java.io.Reader zzReader;

  /** the current state of the DFA */
  private int zzState;

  /** the current lexical state */
  private int zzLexicalState = YYINITIAL;

  /** this buffer contains the current text to be matched and is
      the source of the yytext() string */
  private char zzBuffer[] = new char[ZZ_BUFFERSIZE];

  /** the textposition at the last accepting state */
  private int zzMarkedPos;

  /** the current text position in the buffer */
  private int zzCurrentPos;

  /** startRead marks the beginning of the yytext() string in the buffer */
  private int zzStartRead;

  /** endRead marks the last character in the buffer, that has been read
      from input */
  private int zzEndRead;

  /** number of newlines encountered up to the start of the matched text */
  private int yyline;

  /** the number of characters up to the start of the matched text */
  private int yychar;

  /**
   * the number of characters from the last newline up to the start of the 
   * matched text
   */
  private int yycolumn;

  /** 
   * zzAtBOL == true <=> the scanner is currently at the beginning of a line
   */
  private boolean zzAtBOL = true;

  /** zzAtEOF == true <=> the scanner is at the EOF */
  private boolean zzAtEOF;

  /** denotes if the user-EOF-code has already been executed */
  private boolean zzEOFDone;

  /* user code: */
    static private FrameworkLogger logger =
            FrameworkLogger.getLog(SqlFileScanner.class);
    private StringBuffer commandBuffer = new StringBuffer();
    private boolean interactive;
    private PrintStream psStd = System.out;
    private String magicPrefix;
    private int requestedState = YYINITIAL;
    private String rawLeadinPrompt;
    private boolean specialAppendState;
    // This last is needed for very unique check needed when appending to
    // a SQL command.  Only applies to interactive mode.

    public void setRequestedState(int requestedState) {
        this.requestedState = requestedState;
    }

    /**
     * Really need a way to validate that this is called before using the
     * scanner, like Spring's init-method property.
     * For now, will just check explicitly before using.
     */
    public void setRawLeadinPrompt(String rawLeadinPrompt) {
        this.rawLeadinPrompt = rawLeadinPrompt;
    }

    private void rawLeadinPrompt() {
        if (!interactive) {
            return;
        }
        assert rawLeadinPrompt != null:
            "Internal assertion failed.  "
            + "Scanner's message Resource Bundle not initialized properly";
        psStd.println(rawLeadinPrompt);
    }

    // Trims only the end
    private void trimBuffer() {
        int len = commandBuffer.length();
        commandBuffer.setLength(len -
            ((len > 1 && commandBuffer.charAt(len - 2) == '\r') ? 2 : 1));
    }

    public void setCommandBuffer(String s) {
        commandBuffer.setLength(0);
        commandBuffer.append(s);
    }

    public void setInteractive(boolean interactive) {
        this.interactive = interactive;
    }

    public void setMagicPrefix(String magicPrefix) {
        this.magicPrefix = magicPrefix;
    }

    public void setStdPrintStream(PrintStream psStd) {
        this.psStd = psStd;
    }

    //private String sqlPrompt = "+sql> ";
    private String sqlPrompt = null;
    public void setSqlPrompt(String sqlPrompt)
    {
        this.sqlPrompt = sqlPrompt;
    }
    public String getSqlPrompt() {
        return sqlPrompt;
    }

    //private String sqltoolPrompt = "sql> ";
    private String sqltoolPrompt = null;
    public void setSqltoolPrompt(String sqltoolPrompt)
    {
        this.sqltoolPrompt = sqltoolPrompt;
    }
    public String getSqltoolPrompt() {
        return sqltoolPrompt;
    }
    //private String rawPrompt = "raw> ";
    private String rawPrompt = null;
    public void setRawPrompt(String rawPrompt)
    {
        this.rawPrompt = rawPrompt;
    }
    public String getRawPrompt() {
        return rawPrompt;
    }

    private void debug(String id, String msg) {
        logger.finest(id + ":  [" + msg + ']');
    }

    public String strippedYytext() {
        String lineString = yytext();
        int len = lineString.length();
        len = len - ((len > 1 && lineString.charAt(len - 2) == '\r') ? 2 : 1);
        return (lineString.substring(0, len));
    }

    // Trims only the end
    public void pushbackTrim() {
        String lineString = yytext();
        int len = lineString.length();
        yypushback((len > 1 && lineString.charAt(len - 2) == '\r') ? 2 : 1);
    }

    private void prompt(String s) {
        if (!interactive) return;
        psStd.print(s);
    }

    public void prompt() {
        if (sqltoolPrompt != null) prompt(sqltoolPrompt);
        specialAppendState = (interactive && magicPrefix != null);
        // This tells scanner that if SQL input "looks" empty, it isn't.
        if (interactive && magicPrefix != null) {
            psStd.print(magicPrefix);
            magicPrefix = null;
        }
    }


  /**
   * Creates a new scanner
   * There is also a java.io.InputStream version of this constructor.
   *
   * @param   in  the java.io.Reader to read input from.
   */
  public SqlFileScanner(java.io.Reader in) {
    this.zzReader = in;
  }

  /**
   * Creates a new scanner.
   * There is also java.io.Reader version of this constructor.
   *
   * @param   in  the java.io.Inputstream to read input from.
   */
  public SqlFileScanner(java.io.InputStream in) {
    this(new java.io.InputStreamReader(in));
  }

  /** 
   * Unpacks the compressed character translation table.
   *
   * @param packed   the packed character translation table
   * @return         the unpacked character translation table
   */
  private static char [] zzUnpackCMap(String packed) {
    char [] map = new char[0x10000];
    int i = 0;  /* index in packed string  */
    int j = 0;  /* index in unpacked array */
    while (i < 136) {
      int  count = packed.charAt(i++);
      char value = packed.charAt(i++);
      do map[j++] = value; while (--count > 0);
    }
    return map;
  }


  /**
   * Refills the input buffer.
   *
   * @return      <code>false</code>, iff there was new input.
   * 
   * @exception   java.io.IOException  if any I/O-Error occurs
   */
  private boolean zzRefill() throws java.io.IOException {

    /* first: make room (if you can) */
    if (zzStartRead > 0) {
      System.arraycopy(zzBuffer, zzStartRead,
                       zzBuffer, 0,
                       zzEndRead-zzStartRead);

      /* translate stored positions */
      zzEndRead-= zzStartRead;
      zzCurrentPos-= zzStartRead;
      zzMarkedPos-= zzStartRead;
      zzStartRead = 0;
    }

    /* is the buffer big enough? */
    if (zzCurrentPos >= zzBuffer.length) {
      /* if not: blow it up */
      char newBuffer[] = new char[zzCurrentPos*2];
      System.arraycopy(zzBuffer, 0, newBuffer, 0, zzBuffer.length);
      zzBuffer = newBuffer;
    }

    /* finally: fill the buffer with new input */
    int numRead = zzReader.read(zzBuffer, zzEndRead,
                                            zzBuffer.length-zzEndRead);

    if (numRead > 0) {
      zzEndRead+= numRead;
      return false;
    }
    // unlikely but not impossible: read 0 characters, but not at end of stream    
    if (numRead == 0) {
      int c = zzReader.read();
      if (c == -1) {
        return true;
      } else {
        zzBuffer[zzEndRead++] = (char) c;
        return false;
      }     
    }

	// numRead < 0
    return true;
  }

    
  /**
   * Closes the input stream.
   */
  public final void yyclose() throws java.io.IOException {
    zzAtEOF = true;            /* indicate end of file */
    zzEndRead = zzStartRead;  /* invalidate buffer    */

    if (zzReader != null)
      zzReader.close();
  }


  /**
   * Resets the scanner to read from a new input stream.
   * Does not close the old reader.
   *
   * All internal variables are reset, the old input stream 
   * <b>cannot</b> be reused (internal buffer is discarded and lost).
   * Lexical state is set to <tt>ZZ_INITIAL</tt>.
   *
   * @param reader   the new input stream 
   */
  public final void yyreset(java.io.Reader reader) {
    zzReader = reader;
    zzAtBOL  = true;
    zzAtEOF  = false;
    zzEOFDone = false;
    zzEndRead = zzStartRead = 0;
    zzCurrentPos = zzMarkedPos = 0;
    yyline = yychar = yycolumn = 0;
    zzLexicalState = YYINITIAL;
  }


  /**
   * Returns the current lexical state.
   */
  public final int yystate() {
    return zzLexicalState;
  }


  /**
   * Enters a new lexical state
   *
   * @param newState the new lexical state
   */
  public final void yybegin(int newState) {
    zzLexicalState = newState;
  }


  /**
   * Returns the text matched by the current regular expression.
   */
  public final String yytext() {
    return new String( zzBuffer, zzStartRead, zzMarkedPos-zzStartRead );
  }


  /**
   * Returns the character at position <tt>pos</tt> from the 
   * matched text. 
   * 
   * It is equivalent to yytext().charAt(pos), but faster
   *
   * @param pos the position of the character to fetch. 
   *            A value from 0 to yylength()-1.
   *
   * @return the character at position pos
   */
  public final char yycharat(int pos) {
    return zzBuffer[zzStartRead+pos];
  }


  /**
   * Returns the length of the matched text region.
   */
  public final int yylength() {
    return zzMarkedPos-zzStartRead;
  }


  /**
   * Reports an error that occured while scanning.
   *
   * In a wellformed scanner (no or only correct usage of 
   * yypushback(int) and a match-all fallback rule) this method 
   * will only be called with things that "Can't Possibly Happen".
   * If this method is called, something is seriously wrong
   * (e.g. a JFlex bug producing a faulty scanner etc.).
   *
   * Usual syntax/scanner level error handling should be done
   * in error fallback rules.
   *
   * @param   errorCode  the code of the errormessage to display
   */
  private void zzScanError(int errorCode) {
    String message;
    try {
      message = ZZ_ERROR_MSG[errorCode];
    }
    catch (ArrayIndexOutOfBoundsException e) {
      message = ZZ_ERROR_MSG[ZZ_UNKNOWN_ERROR];
    }

    throw new Error(message);
  } 


  /**
   * Pushes the specified amount of characters back into the input stream.
   *
   * They will be read again by then next call of the scanning method
   *
   * @param number  the number of characters to be read again.
   *                This number must not be greater than yylength()!
   */
  public void yypushback(int number)  {
    if ( number > yylength() )
      zzScanError(ZZ_PUSHBACK_2BIG);

    zzMarkedPos -= number;
  }


  /**
   * Contains user EOF-code, which will be executed exactly once,
   * when the end of file is reached
   */
  private void zzDoEOF() throws java.io.IOException {
    if (!zzEOFDone) {
      zzEOFDone = true;
      yyclose();
    }
  }


  /**
   * Resumes scanning until the next regular expression is matched,
   * the end of input is encountered or an I/O-Error occurs.
   *
   * @return      the next token
   * @exception   java.io.IOException  if any I/O-Error occurs
   */
  public Token yylex() throws java.io.IOException {
    int zzInput;
    int zzAction;

    // cached fields:
    int zzCurrentPosL;
    int zzMarkedPosL;
    int zzEndReadL = zzEndRead;
    char [] zzBufferL = zzBuffer;
    char [] zzCMapL = ZZ_CMAP;

    int [] zzTransL = ZZ_TRANS;
    int [] zzRowMapL = ZZ_ROWMAP;
    int [] zzAttrL = ZZ_ATTRIBUTE;

    while (true) {
      zzMarkedPosL = zzMarkedPos;

      boolean zzR = false;
      for (zzCurrentPosL = zzStartRead; zzCurrentPosL < zzMarkedPosL;
                                                             zzCurrentPosL++) {
        switch (zzBufferL[zzCurrentPosL]) {
        case '\u000B':
        case '\u000C':
        case '\u0085':
        case '\u2028':
        case '\u2029':
          yyline++;
          yycolumn = 0;
          zzR = false;
          break;
        case '\r':
          yyline++;
          yycolumn = 0;
          zzR = true;
          break;
        case '\n':
          if (zzR)
            zzR = false;
          else {
            yyline++;
            yycolumn = 0;
          }
          break;
        default:
          zzR = false;
          yycolumn++;
        }
      }

      if (zzR) {
        // peek one character ahead if it is \n (if we have counted one line too much)
        boolean zzPeek;
        if (zzMarkedPosL < zzEndReadL)
          zzPeek = zzBufferL[zzMarkedPosL] == '\n';
        else if (zzAtEOF)
          zzPeek = false;
        else {
          boolean eof = zzRefill();
          zzEndReadL = zzEndRead;
          zzMarkedPosL = zzMarkedPos;
          zzBufferL = zzBuffer;
          if (eof) 
            zzPeek = false;
          else 
            zzPeek = zzBufferL[zzMarkedPosL] == '\n';
        }
        if (zzPeek) yyline--;
      }
      if (zzMarkedPosL > zzStartRead) {
        switch (zzBufferL[zzMarkedPosL-1]) {
        case '\n':
        case '\u000B':
        case '\u000C':
        case '\u0085':
        case '\u2028':
        case '\u2029':
          zzAtBOL = true;
          break;
        case '\r': 
          if (zzMarkedPosL < zzEndReadL)
            zzAtBOL = zzBufferL[zzMarkedPosL] != '\n';
          else if (zzAtEOF)
            zzAtBOL = false;
          else {
            boolean eof = zzRefill();
            zzMarkedPosL = zzMarkedPos;
            zzEndReadL = zzEndRead;
            zzBufferL = zzBuffer;
            if (eof) 
              zzAtBOL = false;
            else 
              zzAtBOL = zzBufferL[zzMarkedPosL] != '\n';
          }
          break;
        default:
          zzAtBOL = false;
        }
      }
      zzAction = -1;

      zzCurrentPosL = zzCurrentPos = zzStartRead = zzMarkedPosL;
  
      if (zzAtBOL)
        zzState = ZZ_LEXSTATE[zzLexicalState+1];
      else
        zzState = ZZ_LEXSTATE[zzLexicalState];


      zzForAction: {
        while (true) {
    
          if (zzCurrentPosL < zzEndReadL)
            zzInput = zzBufferL[zzCurrentPosL++];
          else if (zzAtEOF) {
            zzInput = YYEOF;
            break zzForAction;
          }
          else {
            // store back cached positions
            zzCurrentPos  = zzCurrentPosL;
            zzMarkedPos   = zzMarkedPosL;
            boolean eof = zzRefill();
            // get translated positions and possibly new buffer
            zzCurrentPosL  = zzCurrentPos;
            zzMarkedPosL   = zzMarkedPos;
            zzBufferL      = zzBuffer;
            zzEndReadL     = zzEndRead;
            if (eof) {
              zzInput = YYEOF;
              break zzForAction;
            }
            else {
              zzInput = zzBufferL[zzCurrentPosL++];
            }
          }
          int zzNext = zzTransL[ zzRowMapL[zzState] + zzCMapL[zzInput] ];
          if (zzNext == -1) break zzForAction;
          zzState = zzNext;

          int zzAttributes = zzAttrL[zzState];
          if ( (zzAttributes & 1) == 1 ) {
            zzAction = zzState;
            zzMarkedPosL = zzCurrentPosL;
            if ( (zzAttributes & 8) == 8 ) break zzForAction;
          }

        }
      }

      // store back cached position
      zzMarkedPos = zzMarkedPosL;

      switch (zzAction < 0 ? zzAction : ZZ_ACTION[zzAction]) {
        case 19: 
          { commandBuffer.append(yytext());
        debug("SQL '", yytext());
        yybegin(SQL);
          }
        case 35: break;
        case 9: 
          { commandBuffer.setLength(0);
    yybegin(SPECIAL);
          }
        case 36: break;
        case 30: 
          { pushbackTrim();
        /* embedded comment may disable opening quotes and closing ; */
        debug("Spl. -- Comment", yytext());
          }
        case 37: break;
        case 10: 
          { commandBuffer.setLength(0);
    yybegin(EDIT);
          }
        case 38: break;
        case 21: 
          { yybegin(YYINITIAL);
    debug("Gobbled", yytext());
    prompt();
          }
        case 39: break;
        case 31: 
          { /* Ignore top-level traditional comments */
    debug ("/**/ Comment", yytext());
          }
        case 40: break;
        case 8: 
          { return new Token(Token.SQL_TYPE, yyline);
          }
        case 41: break;
        case 2: 
          { prompt();
          }
        case 42: break;
        case 34: 
          { /* These are commands which may contain nested commands and/or which
     * require the closing semicolon to send to the DB engine.
     * The BEGIN and DECLARE needed for PL/SQL probably do not need to
     * terminate the line, as we have it specified here, but I'd rather not be
     * too liberal with proprietary SQL like this, because it's easy to
     * envision other proprietary or non-proprietary commands beginning with
     * DECLARE or BEGIN. */
    setCommandBuffer(strippedYytext());
    yybegin(RAW);
    rawLeadinPrompt();
    if (rawPrompt != null) prompt(rawPrompt);
          }
        case 43: break;
        case 22: 
          { if (commandBuffer.toString().trim().equals(".")) {
        commandBuffer.setLength(0);
        yybegin(RAW);
        rawLeadinPrompt();
        if (rawPrompt != null) prompt(rawPrompt);
    } else {
        requestedState = YYINITIAL;
        yybegin(PROMPT_CHANGE_STATE);
        pushbackTrim();
        return new Token(Token.SPECIAL_TYPE, commandBuffer, yyline);
    }
          }
        case 44: break;
        case 28: 
          { specialAppendState = false;
        commandBuffer.append(yytext());
        /* embedded comment may disable opening quotes and closing ; */
        debug("SQL -- Comment", yytext());
          }
        case 45: break;
        case 17: 
          { if (commandBuffer.length() > 0) commandBuffer.append('\n');
        commandBuffer.append(strippedYytext());
        if (rawPrompt != null) prompt(rawPrompt);
          }
        case 46: break;
        case 26: 
          { yybegin(requestedState);
    prompt();
          }
        case 47: break;
        case 4: 
          { commandBuffer.setLength(0);
    yybegin(MACRO);
          }
        case 48: break;
        case 18: 
          { commandBuffer.append(yytext());
          }
        case 49: break;
        case 11: 
          { specialAppendState = false;
        commandBuffer.append(yytext());
          }
        case 50: break;
        case 25: 
          { requestedState = YYINITIAL;
    yybegin(PROMPT_CHANGE_STATE);
    pushbackTrim();
    return new Token(Token.MACRO_TYPE, commandBuffer, yyline);
          }
        case 51: break;
        case 16: 
          { if (interactive && !specialAppendState) {
            requestedState = YYINITIAL;
            yybegin(PROMPT_CHANGE_STATE);
            pushbackTrim();
            trimBuffer();
            return new Token(Token.BUFFER_TYPE, commandBuffer, yyline);
        }
        specialAppendState = false;
        commandBuffer.append(yytext());
          }
        case 52: break;
        case 29: 
          { yybegin(YYINITIAL);
        prompt();
        return new Token(Token.RAWEXEC_TYPE, commandBuffer, yyline);
          }
        case 53: break;
        case 27: 
          { yybegin(YYINITIAL);
        prompt();
        return new Token(Token.RAW_TYPE, commandBuffer, yyline);
          }
        case 54: break;
        case 14: 
          { specialAppendState = false;
        yybegin(YYINITIAL);
        return new Token(Token.SQL_TYPE, commandBuffer, yyline);
          }
        case 55: break;
        case 33: 
          { /* embedded comment may disable opening closing \n */
        debug("Spl. /**/ Comment", yytext());
          }
        case 56: break;
        case 3: 
          { yybegin(GOBBLE);
    return new Token(Token.SYNTAX_ERR_TYPE, yytext(), yyline);
          }
        case 57: break;
        case 20: 
          { commandBuffer.append(yytext());
        yybegin(SQL);
        debug("SQL \"", yytext());
          }
        case 58: break;
        case 1: 
          { setCommandBuffer(yytext());
    yybegin(SQL);
          }
        case 59: break;
        case 23: 
          { requestedState = YYINITIAL;
    yybegin(PROMPT_CHANGE_STATE);
    pushbackTrim();
    return new Token(Token.PL_TYPE, commandBuffer, yyline);
          }
        case 60: break;
        case 6: 
          { /* Ignore top-level white space */
    debug("Whitespace", yytext());
          }
        case 61: break;
        case 12: 
          { specialAppendState = false;
        commandBuffer.append(yytext());
        if (sqlPrompt != null) prompt(sqlPrompt);
          }
        case 62: break;
        case 24: 
          { requestedState = YYINITIAL;
    yybegin(PROMPT_CHANGE_STATE);
    pushbackTrim();
    return new Token(Token.EDIT_TYPE, commandBuffer, yyline);
          }
        case 63: break;
        case 7: 
          { debug ("-- Comment", yytext());
          }
        case 64: break;
        case 15: 
          { specialAppendState = false;
        commandBuffer.append(yytext());
        yybegin(SQL_SINGLE_QUOTED);
          }
        case 65: break;
        case 5: 
          { commandBuffer.setLength(0);
    yybegin(PL);
          }
        case 66: break;
        case 32: 
          { specialAppendState = false;
        commandBuffer.append(yytext());
        /* embedded comment may disable opening quotes and closing ; */
        debug("SQL /**/ Comment", yytext());
          }
        case 67: break;
        case 13: 
          { specialAppendState = false;
        commandBuffer.append(yytext());
        yybegin(SQL_DOUBLE_QUOTED);
          }
        case 68: break;
        default: 
          if (zzInput == YYEOF && zzStartRead == zzCurrentPos) {
            zzAtEOF = true;
            zzDoEOF();
            switch (zzLexicalState) {
            case SPECIAL: {
              yybegin(YYINITIAL);
    return new Token(Token.UNTERM_TYPE, commandBuffer, yyline);
            }
            case 130: break;
            case SQL_DOUBLE_QUOTED: {
              yybegin(YYINITIAL);
    return new Token(Token.UNTERM_TYPE, commandBuffer, yyline);
            }
            case 131: break;
            case SQL_SINGLE_QUOTED: {
              yybegin(YYINITIAL);
    return new Token(Token.UNTERM_TYPE, commandBuffer, yyline);
            }
            case 132: break;
            case RAW: {
              yybegin(YYINITIAL);
    return new Token(Token.UNTERM_TYPE, commandBuffer, yyline);
            }
            case 133: break;
            case SQL: {
              yybegin(YYINITIAL);
    return new Token(Token.UNTERM_TYPE, commandBuffer, yyline);
            }
            case 134: break;
            case EDIT: {
              yybegin(YYINITIAL);
    return new Token(Token.UNTERM_TYPE, commandBuffer, yyline);
            }
            case 135: break;
            case PL: {
              yybegin(YYINITIAL);
    return new Token(Token.UNTERM_TYPE, commandBuffer, yyline);
            }
            case 136: break;
            case MACRO: {
              yybegin(YYINITIAL);
    return new Token(Token.UNTERM_TYPE, commandBuffer, yyline);
            }
            case 137: break;
            default:
            return null;
            }
          } 
          else {
            zzScanError(ZZ_NO_MATCH);
          }
      }
    }
  }


}
