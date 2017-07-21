package org.comdnmr.cpmgfit2.gui;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import javafx.application.Application;
import javafx.scene.Scene;
import javafx.scene.layout.StackPane;
import javafx.stage.Stage;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;

public class MainApp extends Application {

    public static Stage primaryStage;
    public static ScriptEngine engine;

    public static void main(String[] args) {
        engine = new ScriptEngineManager().getEngineByName("python");
        engine.put(ScriptEngine.ARGV, args);
        launch(args);
    }

    @Override
    public void start(Stage primaryStage) {
        this.primaryStage = primaryStage;
        primaryStage.setScene(new Scene(new StackPane(), 200, 200));
        primaryStage.setTitle("CoMD/NMR CPMG Analysis CPMGFit2");
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
//        try (FileReader fileReader = new FileReader(parameters.getRaw().get(0))) {
//            engine.eval(fileReader);
//        } catch (FileNotFoundException fnfE) {
//        } catch (IOException | ScriptException ioE) {
//        }
    }
}
