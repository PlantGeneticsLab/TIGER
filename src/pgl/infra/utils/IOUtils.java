/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pgl.infra.utils;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStreamWriter;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 *
 * @author Fei Lu
 */
public class IOUtils {
    
    public static BufferedReader getTextGzipReader (String infileS) {
        BufferedReader br = null;
        try {
            //br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(infileS))));
            br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(infileS), 65536)), 65536);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return br;
    }
    
    public static BufferedReader getTextGzipReader (String infileS, int bufferSize) {
        BufferedReader br = null;
        try {
            //br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(infileS))));
            br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(infileS), bufferSize)), bufferSize);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return br;
    }
    
    public static BufferedWriter getTextGzipWriter (String outfileS) {
        BufferedWriter bw = null;
        try {
            bw = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfileS), 65536)), 65536);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return bw;
    }
    
    public static BufferedWriter getTextGzipWriter (String outfileS, int bufferSize) {
        BufferedWriter bw = null;
        try {
            bw = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outfileS), bufferSize)), bufferSize);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return bw;
    }
    
    public static BufferedWriter getTextWriter (String outfileS) {
        BufferedWriter bw = null;
        try {
             bw = new BufferedWriter (new FileWriter(outfileS), 65536);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return bw;
    }

    public static BufferedReader getTextReader (String infileS) {
        BufferedReader br = null;
        try {
            br = new BufferedReader (new FileReader(infileS), 65536);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return br;
    }

    public static DataOutputStream getBigBufferedBinaryGzipWriter (String outfileS) {
        DataOutputStream dos = null;
        try {
            dos = new DataOutputStream(new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(outfileS), 1 << 20), 1 << 20));

        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return dos;
    }

    public static DataOutputStream getBinaryGzipWriter (String outfileS) {
        DataOutputStream dos = null;
        try {
            dos = new DataOutputStream(new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(outfileS), 65536), 65536));

        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return dos;
    }

    public static DataOutputStream getBinaryWriter (String outfileS) {
        DataOutputStream dos = null;
        try {
            dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return dos;
    }
    
    public static DataOutputStream getBinaryWriter (String outfileS, int bufferSize) {
        DataOutputStream dos = null;
        try {
            dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), bufferSize));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return dos;
    }

    public static DataInputStream getBinaryGzipReader (String infileS) {
        DataInputStream dis = null;
        try {
            dis = new DataInputStream(new BufferedInputStream(new GZIPInputStream(new FileInputStream(infileS), 65536), 65536));

        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return dis;
    }

    public static DataInputStream getBinaryReader (String infileS) {
        DataInputStream dis = null;
        try {
            dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infileS), 65536));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return dis;
    }
    
    public static DataInputStream getBinaryReader (String infileS, int bufferSize) {
        DataInputStream dis = null;
        try {
            dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infileS), bufferSize));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return dis;
    }
    
    public static ObjectOutputStream getObjectWriter (String outfileS) {
        ObjectOutputStream oos = null;
        try {
            oos = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return oos;
    }
    
    public static ObjectInputStream getObjectReader (String infileS) {
        ObjectInputStream ois = null;
        try {
            ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(infileS), 65536));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return ois;
    }
    
    public static Dyad<FileChannel, ByteBuffer> getNIOChannelBufferReader (String fileS, int bufferSize) {
        FileChannel fc = null;
        ByteBuffer bb = ByteBuffer.allocate(bufferSize);
        try {
            fc = FileChannel.open(Paths.get(fileS), StandardOpenOption.READ);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
//        //another way to get a channel to read
//        try {
//            FileInputStream fis = new FileInputStream(fileS);
//            FileChannel afc = fis.getChannel();
//            ByteBuffer abb = ByteBuffer.allocate(bufferSize);
//        }
//        catch (Exception e) {
//            e.printStackTrace();
//        }
//byte[] bytes = new byte[1024];
//ByteBuffer byteBufferRead = ByteBuffer.wrap(bytes);
        return new Dyad(fc, bb);
    }
    
    public static Dyad<FileChannel, ByteBuffer> getNIOChannelDirectBufferReader (String fileS, int bufferSize) {
        FileChannel fc = null;
        ByteBuffer bb = ByteBuffer.allocateDirect(bufferSize);
        try {
            fc = FileChannel.open(Paths.get(fileS), StandardOpenOption.READ);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return new Dyad(fc, bb);
    }
    
    public static Dyad<FileChannel, ByteBuffer> getNIOChannelBufferWriter (String fileS, int bufferSize) {
        FileChannel fc = null;
        ByteBuffer bb = ByteBuffer.allocate(bufferSize);        
        try {
            fc = FileChannel.open(Paths.get(fileS), StandardOpenOption.CREATE, StandardOpenOption.WRITE, StandardOpenOption.TRUNCATE_EXISTING);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
//        //another way to get a channel to write
//        try {
//            FileOutputStream fos = new FileOutputStream(fileS);
//            FileChannel afc = fos.getChannel();
//            ByteBuffer abb = ByteBuffer.allocate(bufferSize);
//        }
//        catch (Exception e) {
//            e.printStackTrace();
//        }
//                
        return new Dyad(fc, bb);
    }
    
    public static Dyad<FileChannel, ByteBuffer> getNIOChannelDirectBufferWriter (String fileS, int bufferSize) {
        FileChannel fc = null;
        ByteBuffer bb = ByteBuffer.allocateDirect(bufferSize);        
        try {

            fc = FileChannel.open(Paths.get(fileS), StandardOpenOption.CREATE, StandardOpenOption.WRITE, StandardOpenOption.TRUNCATE_EXISTING);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return new Dyad(fc, bb);
    }

    public static List<File> getVisibleFileListInDir (String inDirS) {
        File[] fs = new File(inDirS).listFiles();
        List<File> fList = new ArrayList<>();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            fList.add(fs[i]);
        }
        Collections.sort(fList);
        return fList;
    }

    public static List<File> getDirListInDir (String inDirS) {
        File[] fs = new File(inDirS).listFiles();
        List<File> fList = new ArrayList<>();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isDirectory()) {
                fList.add(fs[i]);
            }
        }
        Collections.sort(fList);
        return fList;
    }

    public static List<File> getDirListInDirStartsWith (String inDirS, String startStr) {
        List<File> fList = getDirListInDir(inDirS);
        List<File> nfList = new ArrayList<>();
        for (int i = 0; i < fList.size(); i++) {
            if (fList.get(i).getName().startsWith(startStr)) nfList.add(fList.get(i));
        }
        return nfList;
    }

    public static List<File> getDirListInDirEndsWith (String inDirS, String endStr) {
        List<File> fList = getDirListInDir(inDirS);
        List<File> nfList = new ArrayList<>();
        for (int i = 0; i < fList.size(); i++) {
            if (fList.get(i).getName().endsWith(endStr)) nfList.add(fList.get(i));
        }
        return nfList;
    }

    public static List<File> getFileListInDir (String inDirS) {
        File[] fs = new File(inDirS).listFiles();
        List<File> fList = new ArrayList<>(Arrays.asList(fs));
        Collections.sort(fList);
        return fList;
    }
    
    public static List<File> getFileListInDirContains (String inDirS, String containStr) {
        File[] fs = new File(inDirS).listFiles();
        fs = listFilesContains(fs, containStr);
        List<File> fList = new ArrayList<>(Arrays.asList(fs));
        Collections.sort(fList);
        return fList;
    }
    
    public static List<File> getFileListInDirStartsWith (String inDirS, String startStr) {
        File[] fs = new File(inDirS).listFiles();
        fs = listFilesStartsWith(fs, startStr);
        List<File> fList = new ArrayList<>(Arrays.asList(fs));
        Collections.sort(fList);
        return fList;
    }
    
    public static List<File> getFileListInDirEndsWith (String inDirS, String endStr) {
        File[] fs = new File(inDirS).listFiles();
        fs = listFilesEndsWith(fs, endStr);
        List<File> fList = new ArrayList<>(Arrays.asList(fs));
        Collections.sort(fList);
        return fList;
    }
    
    public static File[] listFilesContains (File[] fAll, String containStr) {
        List<File> al = new ArrayList();
        for (int i = 0; i < fAll.length; i++) {
            if (fAll[i].getName().contains(containStr)) al.add(fAll[i]);
        }
        return al.toArray(new File[al.size()]);
    }
    
    public static File[] listFilesStartsWith (File[] fAll, String startStr) {
        List<File> al = new ArrayList();
        for (int i = 0; i < fAll.length; i++) {
            if (fAll[i].getName().startsWith(startStr)) al.add(fAll[i]);
        }
        return al.toArray(new File[al.size()]);
    }
    
    public static File[] listFilesEndsWith (File[] fAll, String endStr) {
        List<File> al = new ArrayList();
        for (int i = 0; i < fAll.length; i++) {
            if (fAll[i].getName().endsWith(endStr)) al.add(fAll[i]);
        }
        return al.toArray(new File[al.size()]);
    }
    
    /**
     * List all the files in a directory
     * @param dir
     * @return 
     */
    public static File[] listRecursiveFiles (File dir) {
        TreeSet<File> fSet = getRecursiveFiles (dir);
        return fSet.toArray(new File[fSet.size()]);
    }
    
    private static TreeSet<File> getRecursiveFiles (File dir) {
        TreeSet<File> fileTree = new TreeSet();
        for (File entry : dir.listFiles()) {
            if (entry.isFile()) fileTree.add(entry);
            else fileTree.addAll(getRecursiveFiles(entry));
        }
        return fileTree;
    }
}
