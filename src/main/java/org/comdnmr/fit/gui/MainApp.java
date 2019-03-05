package org.comdnmr.fit.gui;

import java.io.InputStream;
import java.io.InputStreamReader;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.scene.Scene;
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
            Platform.exit();
            System.exit(0);
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
}
