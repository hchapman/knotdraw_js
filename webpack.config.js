module.exports = {
    module: {
        rules: [
            {
                test: /\.js$/,
                exclude: /(js|node_modules|bower_components)/,
                use: {
                    loader: 'babel-loader',
                    options: {
                        presets: ['env']
                    }
                }
            }
        ]
    },
    devtool: "source-map",
    entry: {
        "orthogonal": "./src/orthogonal.js",
        "orth_worker": "./src/orth_worker.js"
    },
    output: {
        filename: "js/[name].js"
    }
}
