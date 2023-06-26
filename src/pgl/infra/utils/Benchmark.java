/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.utils;

/**
 * Measure performance of programs
 * @author Fei Lu
 */
public class Benchmark {
    
    public static long getTimeSpanNanoseconds (long timeStart) {
        return System.nanoTime()-timeStart;
    }
    
    public static double getTimeSpanMilliseconds (long timeStart) {
        return (double)getTimeSpanNanoseconds(timeStart)/1000000;
    }
    
    public static double getTimeSpanSeconds (long timeStart) {
        return (double)getTimeSpanNanoseconds(timeStart)/1000000000;
    }
    
    public static double getTimeSpanMinutes (long timeStart) {
        return (double)getTimeSpanSeconds(timeStart)/60;
    }
    
    public static double getTimeSpanHours (long timeStart) {
        return (double)getTimeSpanMinutes(timeStart)/60;
    }
    
    public static double getUsedMemoryGb () {
        return (double)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024/1024;
    }
}
