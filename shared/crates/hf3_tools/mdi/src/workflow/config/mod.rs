//! mdi::workflow::config:Config stores configuration values in a hash map
//! that can be easily passed to data processing functions as a single variable.

// dependencies
use std::collections::HashMap;
use std::env;
use std::str::FromStr;

/// Macro to define one or more environment variable key constants as `const KEY: &str = "KEY";`.
/// Doing so at the top of a binary or library module is optional but improves code 
/// readability and helps avoid typos in string literals used as keys to access 
/// configuration values, since calls can now take the form `cfg.u8.get(KEY)`, etc.
/// Variable keys set in this way must be all uppercase to follow Rust's naming conventions.
/// Provide keys as either space-separated or comma-separated lists.
#[macro_export]
macro_rules! config_key_constants {
    ($($key:ident)+) => { // support space-separated keys
        $(
            const $key: &str = stringify!($key);
        )+
    };
    ($($key:ident),+ $(,)?) => { // support comma-separated keys    
        $(
            const $key: &str = stringify!($key);
        )+
    };
}

/// The Config struct gathers configuration values from environment variables
/// and can store derived configuration values in a hash map organized by data type.
/// Supported data types include u8 (e.g., boolean 0|1 flags), usize, i32, f64, and String.
pub struct Config {
    pub u8:     HashMap<String, u8>, // including integer 0|1 boolean flags
    pub usize:  HashMap<String, usize>,
    pub i32:    HashMap<String, i32>,
    pub f64:    HashMap<String, f64>,
    pub string: HashMap<String, String>,
    pub bool:   HashMap<String, bool>,
}
impl Config {

    /// Create a new empty Config instance.
    pub fn new() -> Self {
        Config {
            u8:     HashMap::new(),
            usize:  HashMap::new(),
            i32:    HashMap::new(),
            f64:    HashMap::new(),
            string: HashMap::new(),
            bool:   HashMap::new(),
        }
    }

    /* ------------------------------------------------------------------
    environment variable setters
    ------------------------------------------------------------------ */
    /// Set u8 configuration values from environment variables.
    /// Panic if any of the specified keys are not found or cannot be parsed as u8.
    pub fn set_u8_env(&mut self, keys: &[&str]) {
        for &key in keys {
            let value_str = Self::get_env_string(key);
            self.u8.insert(key.to_string(), Self::parse_env_string(key, &value_str, "u8"));
        }
    }
    /// Set usize configuration values from environment variables.
    /// Panic if any of the specified keys are not found or cannot be parsed as usize.
    pub fn set_usize_env(&mut self, keys: &[&str]) {
        for &key in keys {
            let value_str = Self::get_env_string(key);
            self.usize.insert(key.to_string(), Self::parse_env_string(key, &value_str, "usize"));
        }
    }
    /// Set i32 configuration values from environment variables.
    /// Panic if any of the specified keys are not found or cannot be parsed as i32.
    pub fn set_i32_env(&mut self, keys: &[&str]) {
        for &key in keys {
            let value_str = Self::get_env_string(key);
            self.i32.insert(key.to_string(), Self::parse_env_string(key, &value_str, "i32"));
        }
    }
    /// Set f64 configuration values from environment variables.
    /// Panic if any of the specified keys are not found or cannot be parsed as f64.
    pub fn set_f64_env(&mut self, keys: &[&str]) {
        for &key in keys {
            let value_str = Self::get_env_string(key);
            self.f64.insert(key.to_string(), Self::parse_env_string(key, &value_str, "f64"));
        }
    }
    /// Set String configuration values from environment variables.
    /// Panic if any of the specified keys are not found.
    pub fn set_string_env(&mut self, keys: &[&str]) {
        for &key in keys {
            let value_str = Self::get_env_string(key);
            self.string.insert(key.to_string(), value_str);
        }
    }
    /* ------------------------------------------------------------------
    environment variable helpers
    ------------------------------------------------------------------ */
    // get the initial string representation of an environment variable
    fn get_env_string(key: &str) -> String {
        match env::var_os(key) {
            Some(value) => value.to_string_lossy().to_string(),
            None => panic!("Environment variable {key} is not set."),
        }
    }
    // parse an environment variable string into the desired data type
    fn parse_env_string<T: FromStr>(key: &str, value: &str, data_type: &str) -> T {
        match value.parse::<T>() {
            Ok(parsed_value) => parsed_value,
            Err(_) => panic!("Environment variable {key} string value '{value}' could not be parsed as {data_type}."),
        }
    }
    /* ------------------------------------------------------------------
    derived variable setters
    ------------------------------------------------------------------ */
    /// Set a (derived) u8 configuration value directly.
    /// Any existing value is overridden and returned as an Option.
    pub fn set_u8(&mut self, key: &str, value: u8) -> Option<u8> {
        self.u8.insert(key.to_string(), value)
    }
    /// Set a (derived) usize configuration value directly.
    /// Any existing value is overridden and returned as an Option.
    pub fn set_usize(&mut self, key: &str, value: usize) -> Option<usize> {
        self.usize.insert(key.to_string(), value)
    }
    /// Set a (derived) i32 configuration value directly.
    /// Any existing value is overridden and returned as an Option.
    pub fn set_i32(&mut self, key: &str, value: i32) -> Option<i32> {
        self.i32.insert(key.to_string(), value)
    }
    /// Set a (derived) f64 configuration value directly.
    /// Any existing value is overridden and returned as an Option.
    pub fn set_f64(&mut self, key: &str, value: f64) -> Option<f64> {
        self.f64.insert(key.to_string(), value)
    }
    /// Set a (derived) String configuration value directly.
    /// Any existing value is overridden and returned as an Option.
    pub fn set_string(&mut self, key: &str, value: String) -> Option<String> {
        self.string.insert(key.to_string(), value)
    }
    /// Set a (derived) bool configuration value directly.
    /// Any existing value is overridden and returned as an Option.
    pub fn set_bool(&mut self, key: &str, value: bool) -> Option<bool> {
        self.bool.insert(key.to_string(), value)
    }
    /* ------------------------------------------------------------------
    config variable getters
    ------------------------------------------------------------------ */
    /// Get a u8 configuration value by key. Panic if the key is not found.
    pub fn get_u8(&self, key: &str) -> u8 {
        *self.u8.get(key).unwrap_or_else(|| Self::key_not_found(key, "u8"))
    }
    /// Get a usize configuration value by key. Panic if the key is not found.
    pub fn get_usize(&self, key: &str) -> usize {
        *self.usize.get(key).unwrap_or_else(|| Self::key_not_found(key, "usize"))
    }
    /// Get an i32 configuration value by key. Panic if the key is not found.
    pub fn get_i32(&self, key: &str) -> i32 {
        *self.i32.get(key).unwrap_or_else(|| Self::key_not_found(key, "i32"))
    }
    /// Get an f64 configuration value by key. Panic if the key is not found.
    pub fn get_f64(&self, key: &str) -> f64 {
        *self.f64.get(key).unwrap_or_else(|| Self::key_not_found(key, "f64"))
    }
    /// Get a String configuration value by key. Panic if the key is not found.
    pub fn get_string(&self, key: &str) -> &str {
        self.string.get(key).unwrap_or_else(|| Self::key_not_found(key, "String"))
    }
    /// Get a bool configuration value by key. Panic if the key is not found.
    pub fn get_bool(&self, key: &str) -> bool {
        *self.bool.get(key).unwrap_or_else(|| Self::key_not_found(key, "bool"))
    }
    /* ------------------------------------------------------------------
    config getter helpers
    ------------------------------------------------------------------ */
    fn key_not_found<T>(key: &str, data_type: &str) -> T {
        panic!("Config key {key} not found in {data_type} value map.")
    }
}
