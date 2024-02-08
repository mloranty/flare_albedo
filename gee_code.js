////// THIS CODE CALCULATES THE WEEKLY ALBEDO WITHIN SIBERIAN FIRE SCARS 
////// TO RUN PROPERLY, IT NEEDS TO BE RUN SEPERATELY FOR EACH YEAR
////// WRITTEN BY ELIZABETH WEBB, LAST UPDATED FEB 2024


///// DATA PREPARATION

// import fires (this should be a public asset, but if you can't access it, the fires are also available here:
//https://arcticdata.io/catalog/view/doi%3A10.18739%2FA2GB1XJ4M)

var fires = ee.FeatureCollection("users/ewebb/Siberia_fires");

// water mask
var water_mask = ee.ImageCollection('MODIS/006/MOD44W')
                  .filter(ee.Filter.date('2015-01-01', '2015-05-01'))
                  .select('water_mask')

// select white sky albedo product
  var albedo_orig = ee.ImageCollection("MODIS/006/MCD43A3")
     .select("BRDF_Albedo_Band_Mandatory_Quality_shortwave", "Albedo_WSA_shortwave")

     
// clip to only consider the fires and mask out water
// only consider top quality data
 var arcticalbedo = albedo_orig.map(function(img) { 
        return img.clipToCollection(fires)
                    .updateMask(water_mask.first().eq(0))
                    .updateMask(img.select("BRDF_Albedo_Band_Mandatory_Quality_shortwave").eq(0))});

var modisscale = arcticalbedo.first().projection().nominalScale()
                    
/////////////////////                    
/////// FIRST STEP: CREATE AN IMAGE OF THE AVERAGE ALBEDO IN EACH ONE WEEK PERIOD
/////////////////////

//https://gis.stackexchange.com/questions/358740/how-to-average-daily-data-to-create-a-weekly-estimate-in-google-earth-engine
/// Do this by year 
    var startDate = ee.Date('2021-03-01')
    var endDate = ee.Date('2021-09-30')
    var dayOffsets = ee.List.sequence(0, endDate.difference(startDate, 'days').subtract(1),7)
    
    var weeklyMeans = ee.ImageCollection.fromImages(
                          dayOffsets.map(function(dayOffset) {
                            var start = startDate.advance(dayOffset, 'days')
                            var end = start.advance(1, 'week')
                            var year = start.get('year')
                            var dayOfYear = start.getRelative('day', 'year')
                            var date = ee.Date(start).format('yyyy-MM-dd')
                                return arcticalbedo
                                  .filterDate(start, end)
                                  .select('Albedo_WSA_shortwave')
                                  .mean()
                                  .set('date',date) }))  
                            .filter( // filter out empty images
                              ee.Filter.listContains('system:band_names', 'Albedo_WSA_shortwave'))


// https://gis.stackexchange.com/questions/343696/calculate-mean-evi-for-multiple-polygons-across-an-image-collection-in-google-ea
// Define a function to be mapped over the albedo image collection that calculates
// mean albedo per image for all polygons in the 'fires' featureCollection
// using the 'reduceRegions' method.
var ZonalStats = function(image) {
  var meanalbedo = image.reduceRegions({
                          collection: fires,
                          reducer: ee.Reducer.mean(),
                          scale: modisscale});
  // Return the featureCollection with the albedo mean summary per feature, but
  // first...
  return meanalbedo
    // ...remove any features that have a null value for any property.
    .filter(ee.Filter.notNull(['mean']))
    // ...map over the featureCollection to edit properties of each feature.
    .map(function(feature) {
      // Return the feature, but first...
      return feature
        // ...select only two properties of interest.
        //Also rename the 'mean' property returned
        // be default from the above reduceRegions function to 'albedo'.
       .select(['mean', 'UniqueId'], ['albedo', 'UniqueId'])
        // add a date property of the image being reduced.
        .set({'date': image.get('date')});
        //.format('YYYY-MM-dd')});
            });
          };
          
var final = weeklyMeans.map(ZonalStats)
  // Flatten the collection of featureCollections into a single featureCollection.
  .flatten();

//print(final.first())

Export.table.toDrive({
        collection: final,
        folder: 'Albedo_output',
        description: 'Albedo_2021',
         // Specify what to export, to get cleaner result
        selectors: [
          'albedo', 
          'date', 
          'UniqueId' ]
})
