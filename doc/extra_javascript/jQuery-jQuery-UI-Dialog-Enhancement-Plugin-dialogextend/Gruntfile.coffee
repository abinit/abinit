module.exports = (grunt) ->
  grunt.initConfig
    pkg: grunt.file.readJSON('package.json')
    uglify:
      options:
        banner: '/*! <%= pkg.name %> <%= pkg.version %> <%= grunt.template.today("yyyy-mm-dd") %> */\n'
      build:
        src: 'build/jquery.dialogextend.js',
        dest: 'build/jquery.dialogextend.min.js'
    coffee:
      compile:
        files:
          'build/jquery.dialogextend.js':['src/jquery.dialogextend.coffee','src/modules/*.coffee']
  
  grunt.loadNpmTasks 'grunt-contrib-uglify'
  grunt.loadNpmTasks 'grunt-contrib-coffee'
  grunt.registerTask 'default', ['coffee','uglify']