/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pgl.infra.utils;


import com.google.common.base.Splitter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;

/**
 *
 * @author Fei Lu 
 */
public class PStringUtils {
    
    public static Pattern whitespacePattern = Pattern.compile("\\s+");
    
    /**
     * Return a joined String from selected elements of input string array
     * @param input
     * @param index
     * @param delimiter
     * @return 
     */
    public static String join (String[] input, int index[], String delimiter) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < index.length; i++) {
            sb.append(input[index[i]]).append(delimiter);
        }
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }
    
    /**
     * Return a string of number filled with 0 on the left
     * @param n
     * @param num
     * @return 
     */
    public static String getNDigitNumber (int n, int num) {
        String s = String.valueOf(num);
        int cnt = n-s.length();
        for (int i = 0; i < cnt; i++) {
            s = "0"+s;
        }
        return s;
    }
    
    /**
     * Return a string of multiple lines transformed from an one line string
     * @param lineLength
     * @param inputStr
     * @return 
     */
    public static String getMultiplelineString (int lineLength, String inputStr) {
        int left = inputStr.length()%lineLength;
        int n;
        StringBuilder sb = new StringBuilder();
        if (left == 0) {
            n = inputStr.length()/lineLength;
            for (int i = 0; i < n; i++) {
                sb.append(inputStr.substring(i*lineLength, (i+1)*lineLength)).append("\n");
            }
            sb.deleteCharAt(sb.length()-1);
        }
        else {
            n = inputStr.length()/lineLength + 1;
            for (int i = 0; i < n-1; i++) {
                sb.append(inputStr.substring(i*lineLength, (i+1)*lineLength)).append("\n");
            }
            sb.append(inputStr.substring(inputStr.length()-left, inputStr.length()));
        }
        return sb.toString();
    }
    
    /**
     * Return an array of multiple line from an one line string
     * @param lineLength
     * @param inputStr
     * @return 
     */
    public static String[] getMultilineString (int lineLength, String inputStr) {
        int left = inputStr.length()%lineLength;
        int n;
        String[] sub = null;
        if (left == 0) {
            n = inputStr.length()/lineLength;
            sub = new String[n];
            for (int i = 0; i < n; i++) {
                sub[i] = inputStr.substring(i*lineLength, (i+1)*lineLength);
            }
        }
        else {
            n = inputStr.length()/lineLength + 1;
            sub = new String[n];
            for (int i = 0; i < n-1; i++) {
                sub[i] = inputStr.substring(i*lineLength, (i+1)*lineLength);
            }
            sub[n-1] = inputStr.substring(inputStr.length()-left, inputStr.length());
        }
        return sub;
    }
    
    /**
     * Return a list of split String using Guava splitter, delimiter tab
     * @param line
     * @return 
     */
    public static List<String> fastSplit (String line) {
        List<String> ls = fastSplit(line, "\t");
        return ls;
    }
    
    /**
     * Return a list of split String based on white space of any length
     * @param line
     * @return 
     */
    public static List<String> fastSplitOnWhitespace (String line) {
        List<String> ls = new ArrayList<String>();
        Splitter spl = Splitter.on(whitespacePattern);
        Iterator<String> it = spl.split(line).iterator();
        while (it.hasNext()) {
          ls.add(it.next());
        }
        return ls;
    }
    
    /**
     * Return a list of split String using Guava splitter
     * @param line
     * @param splitS
     * @return 
     */
    public static List<String> fastSplit (String line, String splitS) {
        List<String> ls = new ArrayList<String>();
        Splitter spl = Splitter.on(splitS);
        Iterator<String> it = spl.split(line).iterator();
        while (it.hasNext()) {
          ls.add(it.next());
        }
        return ls;
    }
    
    /**
     * Split string and return an array of string
     * @param line
     * @param splitS
     * @param startIndex
     * @param endIndex
     * @return 
     */
    public static List<String> fastSplit(String line, String splitS, int startIndex, int endIndex) {
        return fastSplit(line.substring(startIndex, endIndex), splitS);
    }
}
