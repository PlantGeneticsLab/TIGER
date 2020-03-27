/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl;

import pgl.app.speedCall.SpeedCall;

/**
 *
 * @author feilu
 */
public class PGLAPPEntrance {

    public PGLAPPEntrance (String[] args) {
        //this.speedCall(args);
    }

    public void speedCall (String[] args) {
        new SpeedCall(args);
    }

    public static void main (String[] args) {
        new PGLAPPEntrance(args);
    }
}
