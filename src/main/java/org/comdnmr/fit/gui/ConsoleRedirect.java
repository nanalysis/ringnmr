/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.gui;

import java.io.IOException;
import java.io.OutputStream;
import javafx.application.Platform;

import javafx.scene.control.TextArea;

/**
 * This class extends from OutputStream to redirect output to a TextArea
 *
 * @author Martha Beckwith
 *
 */
public class ConsoleRedirect extends OutputStream {

    private TextArea textArea;

    public ConsoleRedirect(TextArea textArea) {
        this.textArea = textArea;
    }

    @Override
    public void write(int b) throws IOException {
        // redirects data to the text area
        if (Platform.isFxApplicationThread()) {
            textArea.appendText(String.valueOf((char) b));
            // scrolls the text area to the end of data
            //textArea.positionCaret(String.valueOf((char)b).length());
            textArea.appendText("");
        } else {
            Platform.runLater(() -> {
                textArea.appendText(String.valueOf((char) b));
                // scrolls the text area to the end of data
                //textArea.positionCaret(String.valueOf((char)b).length());
                textArea.appendText("");
            });
        }
    }
}
