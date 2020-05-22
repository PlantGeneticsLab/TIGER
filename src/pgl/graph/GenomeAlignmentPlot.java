/*
* To change this license header, choose License Headers in Project Properties.
* To change this template file, choose Tools | Templates
* and open the template in the editor.
*/

package pgl.graph;

import com.itextpdf.awt.DefaultFontMapper;
import com.itextpdf.awt.PdfGraphics2D;
import com.itextpdf.text.Document;
import com.itextpdf.text.Rectangle;
import com.itextpdf.text.pdf.PdfContentByte;
import com.itextpdf.text.pdf.PdfWriter;
import java.awt.Color;
import java.io.FileOutputStream;
import java.util.ArrayList;

/**
 *
 * @author Fei Lu
 */
public class GenomeAlignmentPlot {
    public int queryLength = 0;
    public int hitLength = 0;
    public int[] qStart = null;
    public int[] qEnd = null;
    public int[] hStart = null;
    public int[] hEnd = null;
    public int[] qBrickStart = null;
    public int[] qBrickEnd = null;
    String queryName;
    String hitName;
    int maxBarLength = 1200;
    int width = 1600;
    int height  = 1200;
    int xMargin = 200;
    int yMargin = 200;
    int interval = 500;
    int barHeight = 14;
    
    public GenomeAlignmentPlot (String queryName, String hitName, int queryLength, int hitLength, ArrayList<Integer> qStartList, ArrayList<Integer> qEndList, ArrayList<Integer> hStartList, ArrayList<Integer> hEndList) {
        this.queryName = queryName;
        this.hitName = hitName;
        this.queryLength = queryLength;
        this.hitLength = hitLength;
        qStart = new int[qStartList.size()];
        for (int i = 0; i < qStart.length; i++) qStart[i] = qStartList.get(i);
        qEnd = new int[qEndList.size()];
        for (int i = 0; i < qEnd.length; i++) qEnd[i] = qEndList.get(i);
        hStart = new int[hStartList.size()];
        for (int i = 0; i < hStart.length; i++) hStart[i] = hStartList.get(i);
        hEnd = new int[hEndList.size()];
        for (int i = 0; i < hEnd.length; i++) hEnd[i] = hEndList.get(i);
    }
    
    private void standerdize () {
        double ratio = 0;
        if (queryLength > hitLength) {
            ratio = (double)maxBarLength/queryLength;
        }
        else {
            ratio = (double)maxBarLength/hitLength;
        }
        queryLength = (int)(queryLength*ratio);
        hitLength = (int)(hitLength*ratio);
        for (int i = 0; i < qStart.length; i++) {
            qStart[i] = (int)(qStart[i] * ratio);
            qEnd[i] = (int)(qEnd[i] * ratio);
            hStart[i] = (int)(hStart[i] * ratio);
            hEnd[i] = (int)(hEnd[i] * ratio);
        }
        if (qBrickStart != null) {
            for (int i = 0; i < qBrickStart.length; i++) {
                qBrickStart[i] = (int)(qBrickStart[i] * ratio);
                qBrickEnd[i] = (int)(qBrickEnd[i] * ratio);
                //if (qBrickStart[i] == qBrickEnd[i]) qBrickEnd[i]++;
            }
        }
    }
    
    public void addQueryBrick (int[] qBrickStart, int[] qBrickEnd) {
        this.qBrickStart = qBrickStart;
        this.qBrickEnd = qBrickEnd;
    }
    
    public void addQueryBrick (ArrayList<Integer> qStartList, ArrayList<Integer> qEndList) {
        qBrickStart = new int[qStartList.size()];
        qBrickEnd = new int[qEndList.size()];
        for (int i = 0; i < qBrickStart.length; i++) {
            qBrickStart[i] = qStartList.get(i);
            qBrickEnd[i] = qEndList.get(i);
        }
    }
    
    public void saveGraph (String outfileS) {
        this.standerdize();
        try {
            Document pd = new Document(new Rectangle(width,height));
            PdfWriter pw;
            pw = PdfWriter.getInstance(pd, new FileOutputStream (outfileS));
            pd.open();
            PdfContentByte canvas = pw.getDirectContent();
            DefaultFontMapper mapper = new DefaultFontMapper();
            PdfGraphics2D g2d = new PdfGraphics2D(canvas, width, height, mapper);
            g2d.setColor(Color.cyan);
            g2d.fillRoundRect(xMargin, yMargin, queryLength, barHeight, 15, 15);
            g2d.fillRoundRect(xMargin, yMargin+barHeight+interval, hitLength, barHeight, 15, 15);
            g2d.setColor(Color.black);
            g2d.drawString(queryName, xMargin-100, yMargin);
            g2d.drawString(hitName, xMargin-100, yMargin+barHeight+interval);
            g2d.setColor(Color.blue);
            Color c1 =new Color (255,0,0,150);
            Color c2 =new Color (0,0,255,150);
            for (int i = 0; i < this.hStart.length; i++) {
                int[] xPoints = new int[4];
                int[] yPoints = new int[4];
                xPoints[0] = this.qStart[i]+xMargin;
                yPoints[0] = yMargin+barHeight;
                xPoints[1] = this.qEnd[i]+xMargin;
                yPoints[1] = yMargin+barHeight;
                xPoints[3] = this.hStart[i]+xMargin;
                yPoints[3] = yMargin+barHeight+interval;
                xPoints[2] = this.hEnd[i]+xMargin;
                yPoints[2] = yMargin+barHeight+interval;
                if (hStart[i] < hEnd[i]) g2d.setColor(c1);
                else g2d.setColor(c2);
                g2d.fillPolygon(xPoints, yPoints, 4);
            }
            
            if (qBrickStart != null) {
                g2d.setColor(Color.black);
                for (int i = 0; i < this.qBrickStart.length; i++) {
                    g2d.fillRect(xMargin+qBrickStart[i], yMargin, qBrickEnd[i] - qBrickStart[i]+1, barHeight);
                }
            }
            
            g2d.dispose();
            pd.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
