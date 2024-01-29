/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pgl.graph.r;

/**
 * The interface of graphics using R language.
 * @author Fei Lu
 */
public interface RgraphicsInterface {
    
    /**
     * Set R path in current system.
     * @param rPath 
     */
    public void setRPath (String rPath);
    
    /**
     * Set graph color based R color table，see http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf.
     * @param color 
     */
    public void setColor (String color);
    
    /**
     * Set title of the graph.
     * @param title 
     */
    public void setTitle (String title);
    
    /**
     * Set label of x axis.
     * @param xLab 
     */
    public void setXLabel (String xLab);
    
    /**
     * Set label of y axis.
     * @param yLab 
     */
    public void setYLabel (String yLab);
    
    /**
     * Set title size, default = 1, larger number means larger size, try 1.5, 2.0.
     * @param size 
     */
    public void setTitleSize (double size);
    
    /**
     * Set y axis limit.
     * @param low
     * @param high 
     */
    public void setYLim (double low, double high) ;
    
    /**
     * Set x axis limit.
     * @param low
     * @param high 
     */
    public void setXLim (double low, double high) ;
    
    /**
     * Set position of figure legend.
     * @param x
     * @param y 
     */
    public void setLegendPosition (double x, double y) ;
    
    /**
     * Set the font of title.
     * @param n 1: ordinary, 2: bold, 3: italics
     */
    public void setTitleFont (int n) ;
    
    /**
     * Set XY labels size, default = 1, larger number means larger size, try 1.5, 2.0.
     * @param labSize 
     */
    public void setLableSize (double labSize);
    
    /**
     * Set XY axis size, default = 1, larger number means larger size, try 1.5, 2.0.
     * @param axisSize 
     */
    public void setAxisSize (double axisSize);
    
    /**
     * Set plotting character, "pch", from 0 to 25, default = 1, round circle.
     * @param pch 
     */
    public void setPlottingCharacter (int pch);
    /**
     * Set paper mode, title, labels and axis are smaller.
     */
    public void setPaperMode ();
    
    /**
     * Set slide mode, title, labels and axis are bigger.
     */
    public void setSlideMode ();
    
    /**
     * Set title, labels and axis are bigger.
     */
    public void setLargerSize ();
    
    /**
     * Show graph in a Java panel.
     */
    public void showGraph ();
    
    /**
     * Save graph to pdf file.
     * @param outfileS 
     */
    public void saveGraph (String outfileS);
    
    /**
     * Return the R statement for making plot
     * @return
     */
    public String getPlotStatement ();
}
