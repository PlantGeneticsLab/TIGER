/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pgl.graph.r;

import java.io.File;
import java.util.ArrayList;
import rcaller.RCaller;
import rcaller.RCode;

/**
 * Plot scatter plot with multiple classes, each class with different color, using ggplot2 package
 * @author Fei Lu 
 */
public class ScatterPlotMultiClass extends RgraphicsAbstract {
    double[] x = null;
    double[] y = null;
    String[] classes = null;
    String[] variableNames = null;
    
    public ScatterPlotMultiClass (double[] x, double[] y, String[] classes) {
        this.x = x;
        this.y = y;
        this.classes = classes;
    }
    
    public ScatterPlotMultiClass (String[] variableNames) {
        this.variableNames = variableNames;
    }
    
    public ScatterPlotMultiClass (double[][] xValue, double[][] yValue, String[] variableNames) {
        this.variableNames = variableNames;
        int cnt = xValue[0].length * xValue.length;
        x = new double[cnt];
        y = new double[cnt];
        classes = new String[cnt];
        cnt = 0;
        for (int i = 0; i < xValue.length; i++) {
            for (int j = 0; j < xValue[i].length; j++) {
                x[cnt] = xValue[i][j];
                y[cnt] = yValue[i][j];
                classes[cnt] = variableNames[i];
                cnt++;
            }
        }
    }

    public ScatterPlotMultiClass (int[][] xValue, double[][] yValue, String[] variableNames) {
        this.variableNames = variableNames;
        int cnt = xValue[0].length * xValue.length;
        x = new double[cnt];
        y = new double[cnt];
        classes = new String[cnt];
        cnt = 0;
        for (int i = 0; i < xValue.length; i++) {
            for (int j = 0; j < xValue[i].length; j++) {
                x[cnt] = xValue[i][j];
                y[cnt] = yValue[i][j];
                classes[cnt] = variableNames[i];
                cnt++;
            }
        }
    }
    
    public ScatterPlotMultiClass (ArrayList<double[]> xValueList, ArrayList<double[]> yValueList, String[] variableNames) {
        this.variableNames = variableNames;
        this.setDoubleValues(xValueList, yValueList);
    }
    
    public void setDoubleValues (ArrayList<double[]> xValueList, ArrayList<double[]> yValueList) {
        int cnt = 0;
        for (int i = 0; i < xValueList.size(); i++) cnt+=xValueList.get(i).length;
        x = new double[cnt];
        y = new double[cnt];
        classes = new String[cnt];
        cnt = 0;
        for (int i = 0; i < xValueList.size(); i++) {
            for (int j = 0; j < xValueList.get(i).length; j++) {
                x[cnt] = xValueList.get(i)[j];
                y[cnt] = yValueList.get(i)[j];
                classes[cnt] = variableNames[i];
                cnt++;
            }
        }
    }
    
    public void setIntValues (ArrayList<int[]> xValueList, ArrayList<int[]> yValueList) {
        int cnt = 0;
        for (int i = 0; i < xValueList.size(); i++) cnt+=xValueList.get(i).length;
        x = new double[cnt];
        y = new double[cnt];
        classes = new String[cnt];
        cnt = 0;
        for (int i = 0; i < xValueList.size(); i++) {
            for (int j = 0; j < xValueList.get(i).length; j++) {
                x[cnt] = xValueList.get(i)[j];
                y[cnt] = yValueList.get(i)[j];
                classes[cnt] = variableNames[i];
                cnt++;
            }
        }
    }
    
    @Override
    public void showGraph() {
        RCaller caller = new RCaller();
        caller.setRscriptExecutable(this.rPath);
        RCode rCode = new RCode();
        rCode.addRCode("library(ggplot2)");
        rCode.addDoubleArray("x", x);
        rCode.addDoubleArray("y", y);
        rCode.addStringArray("name", classes);
        String dfStr = "df <- data.frame(Class = factor(name), x = x, y = y)";
        rCode.addRCode(dfStr);
        try {
            File plotFile = rCode.startPlot();
            rCode.addRCode(margin);
            rCode.addRCode(this.getPlotStatement());
            caller.setRCode(rCode);
            caller.runOnly();
            rCode.endPlot();
            rCode.showPlot(plotFile);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Override
    public void saveGraph(String outfileS) {
        outfileS = outfileS.replaceAll("\\\\", "/");
        RCaller caller = new RCaller();
        caller.setRscriptExecutable(this.rPath);
        RCode rCode = new RCode();
        rCode.addRCode("library(ggplot2)");
        rCode.addDoubleArray("x", x);
        rCode.addDoubleArray("y", y);
        rCode.addStringArray("name", classes);
        String dfStr = "df <- data.frame(Class = factor(name), x = x, y = y)";
        rCode.addRCode(dfStr);
        rCode.addRCode("pdf(\""+outfileS+"\")");
        rCode.addRCode(margin);
        rCode.addRCode(this.getPlotStatement());
        rCode.addRCode("dev.off()");
        caller.setRCode(rCode);
        caller.runOnly();
    }

    @Override
    public String getPlotStatement() {
        StringBuilder sb = new StringBuilder ();
        sb.append("qplot(x, y, data = df, geom = \"point\", colour = Class, ");
        sb.append("main=\"").append(title).append("\",");
        sb.append("xlab=\"").append(xLab).append("\",");
        sb.append("ylab=\"").append(yLab).append("\",");
        sb.append("cex.main=").append(titleSize).append(",");
        sb.append("cex.lab=").append(labSize).append(",");
        sb.append("cex.axis=").append(axisSize).append(",");
        sb.append("font.main=").append(titleFont).append(",");
        if (sb.charAt(sb.length()-1) == ',') {
            sb.deleteCharAt(sb.length()-1);
        }
        sb.append(")");
        String statement = sb.toString();
        //System.out.println(statement);
        return statement;
    }
    
}
