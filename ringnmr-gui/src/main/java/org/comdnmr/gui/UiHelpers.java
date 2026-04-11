package org.comdnmr.gui;

import java.io.IOException;
import java.io.InputStream;
import java.io.UncheckedIOException;
import java.net.URI;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.control.Label;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.HBox;
import javafx.util.Duration;

public class UiHelpers {
    private static final String HELPER_MESSAGES_DIR = "helper_messages";

    /** Create a small helper icon that shows a message on hover/focus. */
    public static Label createHelperIcon(String fileName) {
        Label icon = new Label("?");
        icon.getStyleClass().add("helper-icon");

        String message = getMessage(fileName);
        Tooltip tooltip = createTooltip(message);
        Tooltip.install(icon, tooltip);

        // keyboard accessible: show tooltip when focused
        icon.focusedProperty().addListener((obs, oldV, newV) -> {
            if (newV) tooltip.show(
                icon,
                icon
                    .localToScreen(icon.getBoundsInLocal())
                    .getMinX(),
                icon
                    .localToScreen(icon.getBoundsInLocal())
                    .getMaxY()
            );
            else tooltip.hide();
        });

        return icon;
    }

    public static HBox createElementWithHelper(Node element, String fileName) {
        HBox hBox = new HBox();
        hBox.setSpacing(5);
        hBox.setAlignment(Pos.CENTER_LEFT);
        Label helper = createHelperIcon(fileName);
        hBox.getChildren().addAll(element, helper);
        return hBox;
    }

    private static String getMessage(String fileName) {
        String path = "/" + HELPER_MESSAGES_DIR + "/" + fileName;
        try (InputStream in = UiHelpers.class.getResourceAsStream(path)) {
            if (in == null) throw new IllegalArgumentException("Resource not found: " + path);
            return new String(in.readAllBytes(), StandardCharsets.UTF_8);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static Tooltip createTooltip(String message) {
        Tooltip tooltip = new Tooltip(message);
        tooltip.setShowDelay(Duration.ZERO);
        tooltip.setShowDuration(Duration.INDEFINITE);
        tooltip.setHideDelay(Duration.ZERO);
        tooltip.setWrapText(true);
        tooltip.setMaxWidth(300);
        return tooltip;
    }
}
