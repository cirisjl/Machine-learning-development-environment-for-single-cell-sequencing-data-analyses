module.exports = {
    webpack: function (config, env) {
        // Ignore source map warnings from maplibre-gl and plotly
        config.ignoreWarnings = [/Failed to parse source map/];
        return config;
    },
    devServer: function (configFunction) {
        return function (proxy, allowedHost) {
            const config = configFunction(proxy, allowedHost);

            if (config.onBeforeSetupMiddleware || config.onAfterSetupMiddleware) {
                config.setupMiddlewares = (middlewares, devServer) => {
                    if (!devServer) {
                        throw new Error('webpack-dev-server is not defined');
                    }

                    if (config.onBeforeSetupMiddleware) {
                        config.onBeforeSetupMiddleware(devServer);
                    }

                    if (config.onAfterSetupMiddleware) {
                        config.onAfterSetupMiddleware(devServer);
                    }

                    return middlewares;
                };

                delete config.onBeforeSetupMiddleware;
                delete config.onAfterSetupMiddleware;
            }

            return config;
        };
    }
};
