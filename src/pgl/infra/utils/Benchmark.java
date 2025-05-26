/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.utils;

import java.io.FileInputStream;
import java.security.MessageDigest;

/**
 * Measure performance of programs
 * @author Fei Lu
 */
public class Benchmark {

    /**
     * Return timespan of program in nanoseconds
     * @param nanoTimeStart
     * @return
     */
    public static long getTimeSpanNanoseconds (long nanoTimeStart) {
        return System.nanoTime()-nanoTimeStart;
    }

    /**
     * Return timespan of program in milliseconds
     * @param nanoTimeStart
     * @return
     */
    public static double getTimeSpanMilliseconds (long nanoTimeStart) {
        return (double)getTimeSpanNanoseconds(nanoTimeStart)/1000000;
    }

    /**
     * Return timespan of program in seconds
     * @param nanoTimeStart
     * @return
     */
    public static double getTimeSpanSeconds (long nanoTimeStart) {
        return (double)getTimeSpanNanoseconds(nanoTimeStart)/1000000000;
    }

    /**
     * Return timespan of program in minutes
     * @param nanoTimeStart
     * @return
     */
    public static double getTimeSpanMinutes (long nanoTimeStart) {
        return (double)getTimeSpanSeconds(nanoTimeStart)/60;
    }

    /**
     * Return timespan of program in hours
     * @param nanoTimeStart
     * @return
     */
    public static double getTimeSpanHours (long nanoTimeStart) {
        return (double)getTimeSpanMinutes(nanoTimeStart)/60;
    }

    /**
     * Return the used memory size in Gb
     * @return
     */
    public static double getUsedMemoryGb () {
        return (double)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024/1024;
    }

    /**
     * Return the MD5 checksum of a file
     * @param infileS
     * @return
     */
    public static String getMD5Checksum(String infileS) {
        MessageDigest md = null;
        try {
            md = MessageDigest.getInstance("MD5");
            FileInputStream fis = new FileInputStream(infileS);
            byte[] dataBytes = new byte[1024];
            int nread;
            while ((nread = fis.read(dataBytes)) != -1) {
                md.update(dataBytes, 0, nread);
            };
            fis.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        byte[] mdBytes = md.digest();
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < mdBytes.length; i++) {
            sb.append(Integer.toString((mdBytes[i] & 0xff) + 0x100, 16).substring(1));
        }
        return sb.toString();
    }
}
