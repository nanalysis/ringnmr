package org.comdnmr.gui;

import java.util.Optional;

import javafx.scene.control.TextField;
import javafx.scene.input.KeyCode;



@FunctionalInterface
interface InputValidator<T> {
    Optional<T> getValue(String input);
}

abstract class ValidatedTextField<T> extends TextField {

    private ValidationStrategy<T> validator;

    public ValidatedTextField(ValidationStrategy<T> validator) {
        this(validator, "");
    }

    public ValidatedTextField(ValidationStrategy<T> validator, String value) {
        this.validator = validator;
        getStyleClass().add("validated-text-field");

        // Listen for key pressed events
        setOnKeyPressed(event -> {
            if (event.getCode() == KeyCode.ENTER) {
                validateInput();
            }
        });
        setText(value);
        validateInput();
    }

    private void validateInput() {
        boolean valid = getValue().isPresent();
        if (valid) {
            getStyleClass().remove("invalid");
        } else {
            setText("");
            getStyleClass().add("invalid");
        }
    }

    Optional<T> getValue() {
        return validator.getValue(getText());
    }
}
