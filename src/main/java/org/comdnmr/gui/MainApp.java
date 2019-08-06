/*
 * CoMD/NMR Software : A Program for Analyzing NMR Dynamics Data
 * Copyright (C) 2018-2019 Bruce A Johnson
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
package org.comdnmr.gui;

import java.io.InputStream;
import java.io.InputStreamReader;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.scene.Scene;
import javafx.scene.control.Alert;
import javafx.scene.control.ButtonType;
import javafx.scene.layout.StackPane;
import javafx.stage.Stage;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;
import org.python.util.InteractiveInterpreter;

public class MainApp extends Application {

    public static Stage primaryStage;
    public static ScriptEngine engine;
    public static InteractiveInterpreter interpreter = new InteractiveInterpreter();
    public static PreferencesController preferencesController;
    public static ConsoleRedirect console;

    public static void main(String[] args) throws ScriptException {
        engine = new ScriptEngineManager().getEngineByName("jython");
        engine.put(ScriptEngine.ARGV, args);
        engine.eval("from exporter import Writer");
        launch(args);
    }

    @Override
    public void start(Stage primaryStage) {
        this.primaryStage = primaryStage;
        primaryStage.setScene(new Scene(new StackPane(), 200, 200));
        primaryStage.setTitle("CoMD/NMR Dynamics Analysis");
        primaryStage.show();
        Parameters parameters = getParameters();
        try {
            InputStream stream = MainApp.class.getClassLoader().getResourceAsStream("stage.py");
            InputStreamReader streamReader = new InputStreamReader(stream);
            engine.eval(streamReader);
        } catch (ScriptException ioE) {
            System.out.println(ioE.getMessage());
            ioE.printStackTrace();

        }
        primaryStage.setOnCloseRequest(e -> {
            checkExit();
            e.consume();
        });
//        try (FileReader fileReader = new FileReader(parameters.getRaw().get(0))) {
//            engine.eval(fileReader);
//        } catch (FileNotFoundException fnfE) {
//        } catch (IOException | ScriptException ioE) {
//        }
    }

    public static InteractiveInterpreter getInterpreter() {
        return interpreter;
    }

    public static void setConsoleController(ConsoleRedirect controller) {
        console = controller;
    }

    public static void checkExit() {
        Alert alert = new Alert(Alert.AlertType.CONFIRMATION, "Exit program", ButtonType.CANCEL, ButtonType.YES);
        alert.showAndWait().ifPresent(response -> {
            if (response == ButtonType.YES) {
                Platform.exit();
                System.exit(0);
            }
        });

    }
}
