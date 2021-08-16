library(progress)

filling_count = function(person_id, concept_id, all_person, verbose=T){
  stopifnot(length(person_id)==length(concept_id))
  
  uniq_person = unique(all_person)
  uniq_concept = unique(concept_id)
  num_person = length(uniq_person)
  num_concept = length(uniq_concept)
  
  New = matrix(0, nrow = num_person, ncol = num_concept,
               dimnames = list(uniq_person, uniq_concept))
  
  total = ncol(New)
  pb = progress_bar$new(
    format = "[:bar] :current/:total |:percent | :elapsed | ETA::eta",
    total = total, clear = F)
  
  for (i in 1:ncol(New)) {
    if (verbose==T) {pb$tick()}
    feature_person = person_id[concept_id == colnames(New)[i]]
    u_feature_person = names(table(feature_person))
    m = match(u_feature_person, rownames(New))
    New[m, i] = unname(table(feature_person))
  }
  return(New)
}


filling_01_byrow = function(person_id, concept_id, all_person, verbose=T){
  stopifnot(length(person_id)==length(concept_id))
  
  uniq_person = unique(all_person)
  uniq_concept = unique(concept_id)
  num_person = length(uniq_person)
  num_concept = length(uniq_concept)
  
  New = matrix(0, nrow = num_person, ncol = num_concept,
               dimnames = list(uniq_person, uniq_concept))
  
  total = nrow(New)
  pb = progress_bar$new(
    format = "[:bar] :current/:total |:percent | :elapsed | ETA::eta",
    total = total, clear = F)
  
  for(i in 1:nrow(New)){
    if (verbose==T) {pb$tick()}
    person_feature = unique(concept_id[person_id == rownames(New)[i]])
    m = match(person_feature, colnames(New))
    New[i, m] = 1
  }
  return(New)
}


filling_01_bycol = function(person_id, concept_id, all_person, verbose=T){
  stopifnot(length(person_id)==length(concept_id))
  
  uniq_person = unique(all_person)
  uniq_concept = unique(concept_id)
  num_person = length(uniq_person)
  num_concept = length(uniq_concept)
  
  New = matrix(0, nrow = num_person, ncol = num_concept,
               dimnames = list(all_person, uniq_concept))
  
  total = ncol(New)
  pb = progress_bar$new(
    format = "[:bar] :current/:total |:percent | :elapsed | ETA::eta",
    total = total, clear = F)
  
  for(i in 1:ncol(New)){
    if (verbose==T) {pb$tick()}
    feature_person = unique(person_id[concept_id == colnames(New)[i]])
    m = match(feature_person, rownames(New))
    New[m, i] = 1
  }
  return(New)
}