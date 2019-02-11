/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.utils;

/**
 *
 * @author brucejohnson
 */
import java.io.*;
import java.net.*;

public class NMRFxClient {

    final Socket socket;
    final PrintWriter out;
    String hostName = "localhost";

    private NMRFxClient(int port) throws IOException {
        socket = new Socket(hostName, port);
        out = new PrintWriter(socket.getOutputStream(), true);
    }

    public static NMRFxClient makeClient(int port) {
        NMRFxClient client = null;
        try {
            client = new NMRFxClient(port);
        } catch (IOException ioE) {

        }
        return client;
    }

    public void sendMessage(String message) throws IOException {
        out.println(message);
    }
}
