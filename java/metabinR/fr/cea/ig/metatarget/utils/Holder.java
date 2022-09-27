package fr.cea.ig.metatarget.utils;

public class Holder<T> {
    private T value;

    public Holder(T value) {
        setValue(value);
    }

    public T getValue() {
        return value;
    }

    public void setValue(T value) {
        this.value = value;
    }
}
