// Custom JS separate from all the components we're currently using

// Home page search binds to hide presuggestions when overridden by typeahead
$(document).ready(function() {
  $("#home-searchbox-input").bind("focus", function() {
    $(".suggestions").show()
  })
  $("#home-searchbox-input").bind("blur", function() {
    $(".suggestions").hide()
  })
  $("#home-searchbox-input").bind("keyup", function() {
    console.log($(this).val())
    if($(this).val() != "") {
      $(".suggestions").hide()
    } else {
      $(".suggestions").show()
    }
  })
})
